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
// r_light.c

extern "C"
{
#include "../quakedef.h"
}

int	r_dlightframecount;


/*
==================
R_AnimateLight
==================
*/
void R_AnimateLight (void)
{
	int			i,j,k;
	
//
// light animations
// 'm' is normal light, 'a' is no light, 'z' is double bright
	i = (int)(cl.time*10);
	for (j=0 ; j<MAX_LIGHTSTYLES ; j++)
	{
		if (!cl_lightstyle[j].length)
		{
			d_lightstylevalue[j] = 256;
			continue;
		}
		k = i % cl_lightstyle[j].length;
		k = cl_lightstyle[j].map[k] - 'a';
		k = k*22;
		d_lightstylevalue[j] = k;
	}	
}

/*
=============================================================================

DYNAMIC LIGHTS BLEND RENDERING

=============================================================================
*/

static void AddLightBlend (float r, float g, float b, float a2)
{
	float	a;

	v_blend[3] = a = v_blend[3] + a2*(1-v_blend[3]);

	a2 = a2/a;

	v_blend[0] = v_blend[1]*(1-a2) + r*a2;
	v_blend[1] = v_blend[1]*(1-a2) + g*a2;
	v_blend[2] = v_blend[2]*(1-a2) + b*a2;
}

static void R_RenderDlight (dlight_t *light)
{
	int		i, j;
	float	a;
	vec3_t	v;
	float rad;

	rad = light->radius * 0.35;

	VectorSubtract (light->origin, r_origin, v);
	if (Length (v) < rad)
	{	// view is inside the dlight
		AddLightBlend (1, 0.5, 0, light->radius * 0.0003);
		return;
	}

	/*
	glBegin (GL_TRIANGLE_FAN);
	glColor3f (0.2,0.1,0.0);
	for (i=0 ; i<3 ; i++)
		v[i] = light->origin[i] - vpn[i]*rad;
	glVertex3fv (v);
	glColor3f (0,0,0);
	for (i=16 ; i>=0 ; i--)
	{
		a = i/16.0 * M_PI*2;
		for (j=0 ; j<3 ; j++)
			v[j] = light->origin[j] + vright[j]*cosf(a)*rad
				+ vup[j]*sinf(a)*rad;
		glVertex3fv (v);
	}
	glEnd ();
	*/
}

/*
=============
R_RenderDlights
=============
*/
void R_RenderDlights (void)
{
	int		i;
	dlight_t	*l;

	r_dlightframecount = r_framecount + 1;	// because the count hasn't
											//  advanced yet for this frame
	/*glDepthMask (0);
	glDisable (GL_TEXTURE_2D);
	glShadeModel (GL_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc (GL_ONE, GL_ONE);

	l = cl_dlights;
	for (i=0 ; i<MAX_DLIGHTS ; i++, l++)
	{
		if (l->die < cl.time || !l->radius)
			continue;
		R_RenderDlight (l);
	}

	glColor3f (1,1,1);
	glDisable (GL_BLEND);
	glEnable (GL_TEXTURE_2D);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask (1);*/
}


/*
=============================================================================

DYNAMIC LIGHTS

=============================================================================
*/

/*
=============
R_MarkLights
=============
*/
/*
void R_MarkLights (dlight_t *light, int bit, mnode_t *node)
{
	mplane_t	*splitplane;
	float		dist;
	msurface_t	*surf;
	int			i;
	
	if (node->contents < 0)
		return;

	splitplane = node->plane;
	dist = DotProduct (light->origin, splitplane->normal) - splitplane->dist;
	
	if (dist > light->radius)
	{
		R_MarkLights (light, bit, node->children[0]);
		return;
	}
	if (dist < -light->radius)
	{
		R_MarkLights (light, bit, node->children[1]);
		return;
	}
		
// mark the polygons
	surf = cl.worldmodel->surfaces + node->firstsurface;
	for (i=0 ; i<node->numsurfaces ; i++, surf++)
	{
		if (surf->dlightframe != r_dlightframecount)
		{
			surf->dlightbits = 0;
			surf->dlightframe = r_dlightframecount;
		}
		surf->dlightbits |= bit;
	}

	R_MarkLights (light, bit, node->children[0]);
	R_MarkLights (light, bit, node->children[1]);
}
*/

void R_MarkLights (dlight_t *light, int bit, mnode_t *node)
{
	mplane_t	*splitplane;
	float		dist, l, maxdist;
	msurface_t	*surf;
	int		i, j, s, t, sidebit;
	vec3_t		impact;

loc0:
	if (node->contents < 0)
		return;

	splitplane = node->plane;
//	dist = PlaneDiff(light->origin, splitplane);
	dist = DotProduct (light->origin, splitplane->normal) - splitplane->dist;

	if (dist > light->radius)
	{
		node = node->children[0];
		goto loc0;
	}
	if (dist < -light->radius)
	{
		node = node->children[1];
		goto loc0;
	}

	maxdist = light->radius * light->radius;
// mark the polygons
	surf = cl.worldmodel->surfaces + node->firstsurface;
	for (i=0 ; i<node->numsurfaces ; i++, surf++)
	{
		dist = DotProduct (light->origin, surf->plane->normal) - surf->plane->dist;		// JT030305 - fix light bleed through
		if (dist >= 0)
			sidebit = 0;
		else
			sidebit = SURF_PLANEBACK;

		if ( (surf->flags & SURF_PLANEBACK) != sidebit )				//Discoloda
			continue;								//Discoloda

		for (j=0 ; j<3 ; j++)
			impact[j] = light->origin[j] - surf->plane->normal[j]*dist;

		// clamp center of light to corner and check brightness
		l = DotProduct(impact, surf->texinfo->vecs[0]) + surf->texinfo->vecs[0][3] - surf->texturemins[0];
		s = l + 0.5;
		s = bound(0, s, surf->extents[0]);
		s = l - s;
		l = DotProduct(impact, surf->texinfo->vecs[1]) + surf->texinfo->vecs[1][3] - surf->texturemins[1];
		t = l + 0.5;
		t = bound(0, t, surf->extents[1]);
		t = l - t;

		// compare to minimum light
		if ((s*s + t*t + dist*dist) < maxdist)
		{
			if (surf->dlightframe != r_dlightframecount)	// not dynamic until now
			{
				surf->dlightbits = bit;
				surf->dlightframe = r_dlightframecount;
			}
			else	// already dynamic
			{
				surf->dlightbits |= bit;
			}
		}
	}
	if (node->children[0]->contents >= 0)
		R_MarkLights (light, bit, node->children[0]);
	if (node->children[1]->contents >= 0)
		R_MarkLights (light, bit, node->children[1]);
}


/*
=============
R_PushDlights
=============
*/
void R_PushDlights (void)
{
	int		i;
	dlight_t	*l;

	r_dlightframecount = r_framecount + 1;	// because the count hasn't
											//  advanced yet for this frame
	l = cl_dlights;

	for (i=0 ; i<MAX_DLIGHTS ; i++, l++)
	{
		if (l->die < cl.time || !l->radius)
			continue;
		R_MarkLights ( l, 1<<i, cl.worldmodel->nodes );
	}
}


/*
=============================================================================

LIGHT SAMPLING

=============================================================================
*/

mplane_t		*lightplane;
vec3_t			lightspot;

int RecursiveLightPoint (mnode_t *node, vec3_t start, vec3_t end)
{
	int			r;
	float		front, back, frac;
	int			side;
	mplane_t	*plane;
	vec3_t		mid;
	msurface_t	*surf;
	int			s, t, ds, dt;
	int			i;
	mtexinfo_t	*tex;
	byte		*lightmap;
	unsigned	scale;
	int			maps;

	if (node->contents < 0)
		return -1;		// didn't hit anything

// calculate mid point

// FIXME: optimize for axial
	plane = node->plane;
	front = DotProduct (start, plane->normal) - plane->dist;
	back = DotProduct (end, plane->normal) - plane->dist;
	side = front < 0;
	
	if ( (back < 0) == side)
		return RecursiveLightPoint (node->children[side], start, end);

	frac = front / (front-back);
	mid[0] = start[0] + (end[0] - start[0])*frac;
	mid[1] = start[1] + (end[1] - start[1])*frac;
	mid[2] = start[2] + (end[2] - start[2])*frac;
	
// go down front side	
	r = RecursiveLightPoint (node->children[side], start, mid);
	if (r >= 0)
		return r;		// hit something
		
	if ( (back < 0) == side )
		return -1;		// didn't hit anuthing
		
// check for impact on this node
	VectorCopy (mid, lightspot);
	lightplane = plane;

	surf = cl.worldmodel->surfaces + node->firstsurface;
	for (i=0 ; i<node->numsurfaces ; i++, surf++)
	{
		if (surf->flags & SURF_DRAWTILED)
			continue;	// no lightmaps

		tex = surf->texinfo;
		
		s = DotProduct (mid, tex->vecs[0]) + tex->vecs[0][3];
		t = DotProduct (mid, tex->vecs[1]) + tex->vecs[1][3];;

		if (s < surf->texturemins[0] ||
		t < surf->texturemins[1])
			continue;
		
		ds = s - surf->texturemins[0];
		dt = t - surf->texturemins[1];
		
		if ( ds > surf->extents[0] || dt > surf->extents[1] )
			continue;

		if (!surf->samples)
			return 0;

		ds >>= 4;
		dt >>= 4;

		lightmap = surf->samples;
		r = 0;
		if (lightmap)
		{

			lightmap += dt * ((surf->extents[0]>>4)+1) + ds;

			for (maps = 0 ; maps < MAXLIGHTMAPS && surf->styles[maps] != 255 ;
					maps++)
			{
				scale = d_lightstylevalue[surf->styles[maps]];
				r += *lightmap * scale;
				lightmap += ((surf->extents[0]>>4)+1) *
						((surf->extents[1]>>4)+1);
			}
			
			r >>= 8;
		}
		
		return r;
	}

// go down back side
	return RecursiveLightPoint (node->children[!side], mid, end);
}

int R_LightPoint (vec3_t p)
{
	vec3_t		end;
	int			r;
	
	if (!cl.worldmodel->lightdata)
		return 255;
	
	end[0] = p[0];
	end[1] = p[1];
	end[2] = p[2] - 2048;
	
	r = RecursiveLightPoint (cl.worldmodel->nodes, p, end);
	
	if (r == -1)
		r = 0;

	return r;
}


/*
=============================================================================

VERTEX LIGHTING

=============================================================================
*/

/*
float	vlight_pitch = 45;
float	vlight_yaw = 45;
float	vlight_highcut = 128;
float	vlight_lowcut = 60;

#define NUMVERTEXNORMALS	162
extern	float	r_avertexnormals[NUMVERTEXNORMALS][3];

byte	anorm_pitch[NUMVERTEXNORMALS];
byte	anorm_yaw[NUMVERTEXNORMALS];

byte	vlighttable[256][256];

float R_GetVertexLightValue (int index, float apitch, float ayaw)
{
	int	pitchofs, yawofs;
	float	retval;

	pitchofs = anorm_pitch[index] + (apitch * 256 / 360);
	yawofs = anorm_yaw[index] + (ayaw * 256 / 360);
	while (pitchofs > 255)
		pitchofs -= 256;
	while (yawofs > 255)
		yawofs -= 256;
	while (pitchofs < 0)
		pitchofs += 256;
	while (yawofs < 0)
		yawofs += 256;

	retval = vlighttable[pitchofs][yawofs];

	return retval / 256;
}

float R_LerpVertexLight (int index1, int index2, float ilerp, float apitch, float ayaw)
{
	float	lightval1, lightval2, val;

	lightval1 = R_GetVertexLightValue (index1, apitch, ayaw);
	lightval2 = R_GetVertexLightValue (index2, apitch, ayaw);

	val = (lightval2*ilerp) + (lightval1*(1-ilerp));

	return val;
}

void R_ResetAnormTable (void)
{
	int	i, j;
	float	forward, yaw, pitch, angle, sp, sy, cp, cy, precut;
	vec3_t	normal, lightvec;

	// Define the light vector here
	angle = DEG2RAD(vlight_pitch);
	sy = sin(angle);
	cy = cos(angle);
	angle = DEG2RAD(-vlight_yaw);
	sp = sin(angle);
	cp = cos(angle);
	lightvec[0] = cp*cy;
	lightvec[1] = cp*sy;
	lightvec[2] = -sp;

	// First thing that needs to be done is the conversion of the
	// anorm table into a pitch/yaw table

	for (i=0 ; i<NUMVERTEXNORMALS ; i++)
	{
		if (r_avertexnormals[i][1] == 0 && r_avertexnormals[i][0] == 0)
		{
			yaw = 0;
			if (r_avertexnormals[i][2] > 0)
				pitch = 90;
			else
				pitch = 270;
		}
		else
		{
			yaw = (int)(atan2(r_avertexnormals[i][1], r_avertexnormals[i][0]) * 57.295779513082320);
			if (yaw < 0)
				yaw += 360;
	
			forward = sqrt(r_avertexnormals[i][0]*r_avertexnormals[i][0] + r_avertexnormals[i][1]*r_avertexnormals[i][1]);
			pitch = (int)(atan2(r_avertexnormals[i][2], forward) * 57.295779513082320);
			if (pitch < 0)
				pitch += 360;
		}
		anorm_pitch[i] = pitch * 256 / 360;
		anorm_yaw[i] = yaw * 256 / 360;
	}

	// Next, a light value table must be constructed for pitch/yaw offsets
	// DotProduct values

	// DotProduct values never go higher than 2, so store bytes as
	// (product * 127.5)

	for (i=0 ; i<256 ; i++)
	{
		angle = DEG2RAD(i * 360 / 256);
		sy = sin(angle);
		cy = cos(angle);
		for (j=0 ; j<256 ; j++)
		{
			angle = DEG2RAD(j * 360 / 256);
			sp = sin(angle);
			cp = cos(angle);

			normal[0] = cp*cy;
			normal[1] = cp*sy;
			normal[2] = -sp;

			precut = ((DotProduct(normal, lightvec) + 2) * 31.5);
			precut = (precut - (vlight_lowcut)) * 256 / (vlight_highcut - vlight_lowcut);
			if (precut > 255)
				precut = 255;
			if (precut < 0)
				precut = 0;
			vlighttable[i][j] = precut;
		}
	}
}

void R_InitVertexLights (void)
{
	R_ResetAnormTable ();
}

*/
