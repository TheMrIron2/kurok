/*
Copyright (C) 1996-1997 Id Software, Inc.
Copyright (C) 2007 Peter Mackay and Chris Swindle.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/
// gl_warp.c -- sky and water polygons

#include <pspgu.h>

extern "C"
{
#include "../quakedef.h"
}

#include "clipping.hpp"

using namespace quake;

extern	model_t	*loadmodel;

/*int		skytexturenum;*/

int		solidskytexture	= -1;
int		alphaskytexture	= -1;
float	speedscale;		// for top sky and bottom sky

msurface_t	*warpface;

extern cvar_t gl_subdivide_size;

static void BoundPoly (int numverts, float *verts, vec3_t mins, vec3_t maxs)
{
	int		i, j;
	float	*v;

	mins[0] = mins[1] = mins[2] = 9999;
	maxs[0] = maxs[1] = maxs[2] = -9999;
	v = verts;
	for (i=0 ; i<numverts ; i++)
		for (j=0 ; j<3 ; j++, v++)
		{
			if (*v < mins[j])
				mins[j] = *v;
			if (*v > maxs[j])
				maxs[j] = *v;
		}
}

static void SubdividePolygon (int numverts, float *verts)
{
	int		i, j, k;
	vec3_t	mins, maxs;
	float	m;
	float	*v;
	vec3_t	front[64], back[64];
	int		f, b;
	float	dist[64];
	float	frac;
	glpoly_t	*poly;
	float	s, t;

	if (numverts > 60)
		Sys_Error ("numverts = %i", numverts);

	BoundPoly (numverts, verts, mins, maxs);

	for (i=0 ; i<3 ; i++)
	{
		m = (mins[i] + maxs[i]) * 0.5;
		m = gl_subdivide_size.value * floorf (m/gl_subdivide_size.value + 0.5);
		if (maxs[i] - m < 8)
			continue;
		if (m - mins[i] < 8)
			continue;

		// cut it
		v = verts + i;
		for (j=0 ; j<numverts ; j++, v+= 3)
			dist[j] = *v - m;

		// wrap cases
		dist[j] = dist[0];
		v-=i;
		VectorCopy (verts, v);

		f = b = 0;
		v = verts;
		for (j=0 ; j<numverts ; j++, v+= 3)
		{
			if (dist[j] >= 0)
			{
				VectorCopy (v, front[f]);
				f++;
			}
			if (dist[j] <= 0)
			{
				VectorCopy (v, back[b]);
				b++;
			}
			if (dist[j] == 0 || dist[j+1] == 0)
				continue;
			if ( (dist[j] > 0) != (dist[j+1] > 0) )
			{
				// clip point
				frac = dist[j] / (dist[j] - dist[j+1]);
				for (k=0 ; k<3 ; k++)
					front[f][k] = back[b][k] = v[k] + frac*(v[3+k] - v[k]);
				f++;
				b++;
			}
		}

		SubdividePolygon (f, front[0]);
		SubdividePolygon (b, back[0]);
		return;
	}

	poly = static_cast<glpoly_t*>(Hunk_Alloc (sizeof(glpoly_t) + (numverts - 1) * sizeof(glvert_t)));
	poly->next = warpface->polys;
	warpface->polys = poly;
	poly->numverts = numverts;
	for (i=0 ; i<numverts ; i++, verts+= 3)
	{
		VectorCopy (verts, poly->verts[i].xyz);
		s = DotProduct (verts, warpface->texinfo->vecs[0]);
		t = DotProduct (verts, warpface->texinfo->vecs[1]);
		poly->verts[i].st[0] = s;
		poly->verts[i].st[1] = t;
	}
}

/*
================
GL_SubdivideSurface

Breaks a polygon up along axial 64 unit
boundaries so that turbulent and sky warps
can be done reasonably.
================
*/
void GL_SubdivideSurface (msurface_t *fa)
{
	vec3_t		verts[64];
	int			numverts;
	int			i;
	int			lindex;
	float		*vec;

	warpface = fa;

	//
	// convert edges back to a normal polygon
	//
	numverts = 0;
	for (i=0 ; i<fa->numedges ; i++)
	{
		lindex = loadmodel->surfedges[fa->firstedge + i];

		if (lindex > 0)
			vec = loadmodel->vertexes[loadmodel->edges[lindex].v[0]].position;
		else
			vec = loadmodel->vertexes[loadmodel->edges[-lindex].v[1]].position;
		VectorCopy (vec, verts[numverts]);
		numverts++;
	}

	SubdividePolygon (numverts, verts[0]);
}

//=========================================================



// speed up sin calculations - Ed
float	turbsin[] =
{
	#include "../gl_warp_sin.h"
};

/*
=============
EmitWaterPolys

Does a water warp on the pre-fragmented glpoly_t chain
=============
*/
void EmitWaterPolys (msurface_t *fa)
{
	const float real_time	= static_cast<float>(realtime);
	const float scale		= (1.0f / 64);
	const float turbscale	= (256.0f / (2.0f * static_cast<float>(M_PI)));

	// For each polygon...
	for (const glpoly_t* p = fa->polys; p; p = p->next)
	{
		// Allocate memory for this polygon.
		const int		unclipped_vertex_count	= p->numverts;
		glvert_t* const	unclipped_vertices		=
			static_cast<glvert_t*>(sceGuGetMemory(sizeof(glvert_t) * unclipped_vertex_count));

		// Generate each vertex.
		const glvert_t*	src			= p->verts;
		const glvert_t*	last_vertex = src + unclipped_vertex_count;
		glvert_t*		dst			= unclipped_vertices;

		while (src != last_vertex)
		{
			// Get the input UVs.
			const float	os = src->st[0];
			const float	ot = src->st[1];

			// Fill in the vertex data.
			dst->st[0] = (os + turbsin[(int) ((ot * 0.025f + real_time) * turbscale) & 255]) * scale;
			dst->st[1] = (ot + turbsin[(int) ((os * 0.025f + real_time) * turbscale) & 255]) * scale;
			dst->xyz[0] = src->xyz[0];
			dst->xyz[1] = src->xyz[1];
			dst->xyz[2] = src->xyz[2];

			// Next vertex.
			++src;
			++dst;
		}

		// Do these vertices need clipped?
		if (clipping::is_clipping_required(unclipped_vertices, unclipped_vertex_count))
		{
			// Clip the polygon.
			const glvert_t*	clipped_vertices;
			std::size_t		clipped_vertex_count;
			clipping::clip(
				unclipped_vertices,
				unclipped_vertex_count,
				&clipped_vertices,
				&clipped_vertex_count);

			// Any vertices left?
			if (clipped_vertex_count)
			{
				// Copy the vertices to the display list.
				const std::size_t buffer_size = clipped_vertex_count * sizeof(glvert_t);
				glvert_t* const display_list_vertices = static_cast<glvert_t*>(sceGuGetMemory(buffer_size));
				memcpy(display_list_vertices, clipped_vertices, buffer_size);

				// Draw the clipped vertices.
				sceGuDrawArray(
					GU_TRIANGLE_FAN,
					GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
					clipped_vertex_count, 0, display_list_vertices);
			}
		}
		else
		{
			// Draw the vertices.
			sceGuDrawArray(
				GU_TRIANGLE_FAN,
				GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
				unclipped_vertex_count, 0, unclipped_vertices);
		}
	}
}


/*
=============
EmitReflectivePolys

Does a reflective warp on the pre-fragmented glpoly_t chain
=============
*/
void EmitReflectivePolys (msurface_t *fa)
{
	// For each polygon...
	for (const glpoly_t* p = fa->polys; p; p = p->next)
	{
		// Allocate memory for this polygon.
		const int		unclipped_vertex_count	= p->numverts;
		glvert_t* const	unclipped_vertices		=
			static_cast<glvert_t*>(sceGuGetMemory(sizeof(glvert_t) * unclipped_vertex_count));

		// Generate each vertex.
		const glvert_t*	src			= p->verts;
		const glvert_t*	last_vertex = src + unclipped_vertex_count;
		glvert_t*		dst			= unclipped_vertices;

		while (src != last_vertex)
		{
			vec3_t	dir;
			VectorSubtract(src->xyz, r_origin, dir);
			dir[2] *= 3;	// flatten the sphere

			const float length = 6 * 63 / sqrtf(DotProduct(dir, dir));

			dir[0] *= length;
			dir[1] *= length;

			dst->st[0] = (dir[0]) * (1.0f / 256.0f);
			dst->st[1] = (dir[1]) * (1.0f / 256.0f);
			dst->xyz[0] = src->xyz[0];
			dst->xyz[1] = src->xyz[1];
			dst->xyz[2] = src->xyz[2];

			// Next vertex.
			++src;
			++dst;
		}

		// Do these vertices need clipped?
		if (clipping::is_clipping_required(unclipped_vertices, unclipped_vertex_count))
		{
			// Clip the polygon.
			const glvert_t*	clipped_vertices;
			std::size_t		clipped_vertex_count;
			clipping::clip(
				unclipped_vertices,
				unclipped_vertex_count,
				&clipped_vertices,
				&clipped_vertex_count);

			// Any vertices left?
			if (clipped_vertex_count)
			{
				// Copy the vertices to the display list.
				const std::size_t buffer_size = clipped_vertex_count * sizeof(glvert_t);
				glvert_t* const display_list_vertices = static_cast<glvert_t*>(sceGuGetMemory(buffer_size));
				memcpy(display_list_vertices, clipped_vertices, buffer_size);

				// Draw the clipped vertices.
				sceGuDrawArray(
					GU_TRIANGLE_FAN,
					GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
					clipped_vertex_count, 0, display_list_vertices);
			}
		}
		else
		{
			// Draw the vertices.
			sceGuDrawArray(
				GU_TRIANGLE_FAN,
				GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
				unclipped_vertex_count, 0, unclipped_vertices);
		}
	}
}

/*
=============
EmitSkyPolys
=============
*/
void EmitSkyPolys (msurface_t *fa)
{
     
  if (!r_skyclip.value)
  {

	for (const glpoly_t* p = fa->polys; p; p = p->next)
	{
        
		glvert_t* const vertices = static_cast<glvert_t*>(sceGuGetMemory(sizeof(glvert_t) * p->numverts));

		const glvert_t*	src			= p->verts;
		const glvert_t*	last_vertex = src + p->numverts;
		glvert_t*		dst			= vertices;
		
		while (src != last_vertex)
		{
			vec3_t	dir;
			VectorSubtract(src->xyz, r_origin, dir);
			dir[2] *= 3;	// flatten the sphere

			const float length = 6 * 63 / sqrtf(DotProduct(dir, dir));

			dir[0] *= length;
			dir[1] *= length;

			dst->st[0] = (speedscale + dir[0]) * (1.0f / 128.0f);
			dst->st[1] = (speedscale + dir[1]) * (1.0f / 128.0f);
			dst->xyz[0] = src->xyz[0];
			dst->xyz[1] = src->xyz[1];
			dst->xyz[2] = src->xyz[2];

			// Next vertex.
			++src;
			++dst;
		}

		sceGuDrawArray(GU_TRIANGLE_FAN, GU_TEXTURE_32BITF | GU_VERTEX_32BITF, p->numverts, 0, vertices);
	}
  }
  else
  {

    for (const glpoly_t* p = fa->polys; p; p = p->next)
	{
		
		// Allocate memory for this polygon.
		const int		unclipped_vertex_count	= p->numverts;
		glvert_t* const	unclipped_vertices		=
			static_cast<glvert_t*>(sceGuGetMemory(sizeof(glvert_t) * unclipped_vertex_count));

		// Generate each vertex.
		const glvert_t*	src			= p->verts;
		const glvert_t*	last_vertex = src + unclipped_vertex_count;
		glvert_t*		dst			= unclipped_vertices;

		while (src != last_vertex)
		{
			vec3_t	dir;
			VectorSubtract(src->xyz, r_origin, dir);
			dir[2] *= 3;	// flatten the sphere

			const float length = 6 * 63 / sqrtf(DotProduct(dir, dir));

			dir[0] *= length;
			dir[1] *= length;

			dst->st[0] = (speedscale + dir[0]) * (1.0f / 128.0f);
			dst->st[1] = (speedscale + dir[1]) * (1.0f / 128.0f);
			dst->xyz[0] = src->xyz[0];
			dst->xyz[1] = src->xyz[1];
			dst->xyz[2] = src->xyz[2];

			// Next vertex.
			++src;
			++dst;
		}

		// Do these vertices need clipped?
		if (clipping::is_clipping_required(unclipped_vertices, unclipped_vertex_count))
		{
			// Clip the polygon.
			const glvert_t*	clipped_vertices;
			std::size_t		clipped_vertex_count;
			clipping::clip(
				unclipped_vertices,
				unclipped_vertex_count,
				&clipped_vertices,
				&clipped_vertex_count);

			// Any vertices left?
			if (clipped_vertex_count)
			{
				// Copy the vertices to the display list.
				const std::size_t buffer_size = clipped_vertex_count * sizeof(glvert_t);
				glvert_t* const display_list_vertices = static_cast<glvert_t*>(sceGuGetMemory(buffer_size));
				memcpy(display_list_vertices, clipped_vertices, buffer_size);

				// Draw the clipped vertices.
				sceGuDrawArray(
					GU_TRIANGLE_FAN,
					GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
					clipped_vertex_count, 0, display_list_vertices);
			}
		}

		else
		{
			// Draw the vertices.
			sceGuDrawArray(
				GU_TRIANGLE_FAN,
				GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
				unclipped_vertex_count, 0, unclipped_vertices);
        }
	}
    }
}

/*
===============
EmitBothSkyLayers

Does a sky warp on the pre-fragmented glpoly_t chain
This will be called for brushmodels, the world
will have them chained together.
===============
*/
void EmitBothSkyLayers (msurface_t *fa)
{
	GL_Bind (solidskytexture);

	if (kurok)
		speedscale = realtime*2;
	else
		speedscale = realtime*8;

	speedscale -= (int)speedscale & ~127 ;

	EmitSkyPolys (fa);

	sceGuEnable(GU_BLEND);
	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);

	GL_Bind (alphaskytexture);

	if (kurok)
		speedscale = realtime*4;
	else
		speedscale = realtime*16;

	speedscale -= (int)speedscale & ~127 ;

	EmitSkyPolys (fa);

	sceGuDisable(GU_BLEND);
	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGB);
}

/*
=================
R_DrawSkyChain
=================
*/
void R_DrawSkyChain (msurface_t *s)
{
	msurface_t	*fa;

	// used when gl_texsort is on
	GL_Bind(solidskytexture);

	if (kurok)
		speedscale = realtime*2;
	else
		speedscale = realtime*8;

	speedscale -= (int)speedscale & ~127 ;

	for (fa=s ; fa ; fa=fa->texturechain)
		EmitSkyPolys (fa);

	sceGuEnable(GU_BLEND);
	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);

	GL_Bind (alphaskytexture);

	if (kurok)
		speedscale = realtime*4;
	else
		speedscale = realtime*16;

	speedscale -= (int)speedscale & ~127 ;

	for (fa=s ; fa ; fa=fa->texturechain)
		EmitSkyPolys (fa);

	sceGuDisable(GU_BLEND);
	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGB);
}


//===============================================================

/*
=============
R_InitSky

A sky texture is 256*128, with the right side being a masked overlay
==============
*/
//void R_InitSky (texture_t *mt)
void R_InitSky (byte *mt)
{
	byte	trans[128*128];

	const byte* const src = (byte *)mt; //  + mt->offsets[0];

	for (int i=0 ; i<128 ; i++)
	{
		for (int j=0 ; j<128 ; j++)
		{
			const byte p = src[i*256 + j + 128];
			trans[(i*128) + j] = p;
		}
	}

	if (solidskytexture == -1)
		solidskytexture = GL_LoadTexture("", 128, 128, trans, qfalse, GU_LINEAR, 0);

	for (int i=0 ; i<128 ; i++)
	{
		for (int j=0 ; j<128 ; j++)
		{
			const byte p = src[i*256 + j];
			if (p == 0)
				trans[(i*128) + j] = 255;
			else
				trans[(i*128) + j] = p;
		}
	}

	if (alphaskytexture == -1)
		alphaskytexture = GL_LoadTexture("", 128, 128, trans, qfalse, GU_LINEAR, 0);
}

