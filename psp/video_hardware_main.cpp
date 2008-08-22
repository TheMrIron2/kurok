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
// r_main.c

extern "C"
{
#include "../quakedef.h"
#include "video_vertex_lighting.h"
}

#include <pspgu.h>
#include <pspgum.h>

#include "clipping.hpp"

using namespace quake;

entity_t	r_worldentity;

qboolean	r_cache_thrash;		// compatability

vec3_t		modelorg, r_entorigin;
entity_t	*currententity;

int			r_visframecount;	// bumped when going to a new PVS
int			r_framecount;		// used for dlight push checking

mplane_t	frustum[4];

int			c_brush_polys, c_alias_polys;

qboolean	envmap;				// qtrue during envmap command capture

int			currenttexture = -1;		// to avoid unnecessary texture sets

int			cnttextures[2] = {-1, -1};     // cached

int			particletexture;	// little dot for particles
int			playertextures;		// up to 16 color translated skins

int			mirrortexturenum;	// quake texturenum, not gltexturenum
qboolean	mirror;
mplane_t	*mirror_plane;

bool	force_fullbright;
bool    additive;
bool    filter;
bool    alphafunc;
bool    alphafunc2;
bool    fixlight;

//
// view origin
//
vec3_t	vup;
vec3_t	vpn;
vec3_t	vright;
vec3_t	r_origin;

ScePspFMatrix4	r_world_matrix;
ScePspFMatrix4	r_base_world_matrix;

//
// screen size info
//
refdef_t	r_refdef;

mleaf_t		*r_viewleaf, *r_oldviewleaf;

texture_t	*r_notexture_mip;

int		d_lightstylevalue[256];	// 8.8 fraction of base light value

void R_MarkLeaves (void);

cvar_t  r_skyclip = {"r_skyclip", "2560", qtrue}; 				// Adjust skybox clipping distance
cvar_t  r_norefresh = {"r_norefresh","0"};      				// Unknown
cvar_t  r_drawentities = {"r_drawentities","1"}; 				// Toggle entity drawing
cvar_t  r_drawviewmodel = {"r_drawviewmodel","1"};              // Toggle view model drawing
cvar_t  r_speeds = {"r_speeds","0"};                            // Toggle speed statistical display
cvar_t  r_fullbright = {"r_fullbright","0"};                    // Toggle fullbright world and entities
cvar_t  r_lightmap = {"r_lightmap","0"};                        // Toggle drawing of lightmaps only
cvar_t  r_shadows = {"r_shadows","0"};                          // Toggle fake stencil shadows for entities
cvar_t  r_mirroralpha = {"r_mirroralpha","1"};                  // Adjust alpha value of mirror textures (not functional)
cvar_t  r_glassalpha = {"r_glassalpha","0.25"};                 // Adjust alpha value of glass textures
cvar_t  r_wateralpha = {"r_wateralpha","0.6"};                  // Adjust alpha value of liquid textures
cvar_t  r_mipmaps = {"r_mipmaps","0",qtrue};                        // Toggle mip map optimising
cvar_t  r_mipmaps_func = {"r_mipmaps_func","2",qtrue};          // Adjust mip map calculations
cvar_t  r_mipmaps_bias = {"r_mipmaps_bias","-7",qtrue};         // Adjust mip map bias level
cvar_t  r_dynamic = {"r_dynamic","1", qtrue};                           // Toggle lightmap dynamic lighting
cvar_t  r_novis = {"r_novis","0"};                                              // Toggle visibility calculations
cvar_t  r_tex_scale_down = {"r_tex_scale_down","0", qtrue};     // Toggle odd texture size resolution down scaling
cvar_t  r_particles_simple = {"r_particles_simple","0", qtrue}; // Toggle classic particle effect
cvar_t  gl_keeptjunctions = {"gl_keeptjunctions","1", qtrue};          // Keeps edge t-junctions in world aligned.
cvar_t  r_vsync = {"r_vsync","0", qtrue};                       // Toggle V-Sync
cvar_t  r_dithering = {"r_dithering","0", qtrue};               // Toggle hardware color dithering
cvar_t  r_antialias = {"r_antialias","0", qtrue};               // Toggle fake anti alias effect
cvar_t  r_menufade = {"r_menufade","0.5", qtrue}; 				// Adjust menu background alpha
cvar_t  r_showtris = {"r_showtris","0", qfalse}; 					// Display wireframes
cvar_t  r_test = {"r_test","0", qfalse}; 						// developer temp test cvar

cvar_t  r_i_model_animation = { "r_i_model_animation", "1", qtrue}; // Toggle smooth model animation
cvar_t  r_i_model_transform = { "r_i_model_transform", "1", qtrue}; // Toggle smooth model movement
cvar_t  r_model_contrast = { "r_model_contrast", "0", qtrue};       // Toggle high contrast model lighting
cvar_t  r_model_brightness = { "r_model_brightness", "1", qtrue};   // Toggle high brightness model lighting
cvar_t  r_skyfog = {"r_skyfog", "1", qtrue};                        // Toggle sky fog

cvar_t  vlight = {"vl_light", "0", qtrue};
cvar_t  vlight_pitch = {"vl_pitch", "45", qtrue};
cvar_t  vlight_yaw = {"vl_yaw", "45", qtrue};
cvar_t  vlight_highcut = {"vl_highcut", "128", qtrue};
cvar_t  vlight_lowcut = {"vl_lowcut", "60", qtrue};

/*
cvar_t	gl_finish = {"gl_finish","0"};
cvar_t	gl_clear = {"gl_clear","0"};
cvar_t	gl_cull = {"gl_cull","1"};
cvar_t	gl_texsort = {"gl_texsort","1"};
cvar_t	gl_smoothmodels = {"gl_smoothmodels","1"};
cvar_t	gl_affinemodels = {"gl_affinemodels","0"};
cvar_t	gl_polyblend = {"gl_polyblend","1"};
cvar_t	gl_flashblend = {"gl_flashblend","1"};
cvar_t	gl_playermip = {"gl_playermip","0"};
cvar_t	gl_nocolors = {"gl_nocolors","0"};
cvar_t	gl_reporttjunctions = {"gl_reporttjunctions","0"};
cvar_t	gl_doubleeyes = {"gl_doubleeys", "1"};

extern	cvar_t	gl_ztrick;*/

/*
=================
R_CullBox

Returns qtrue if the box is completely outside the frustom
=================
*/
qboolean R_CullBox (vec3_t mins, vec3_t maxs)
{
	int		i;

	for (i=0 ; i<4 ; i++)
		if (BoxOnPlaneSide (mins, maxs, &frustum[i]) == 2)
			return qtrue;
	return qfalse;
}


void R_RotateForEntity (entity_t *e)
{
	// Translate.
	const ScePspFVector3 translation = {
		e->origin[0], e->origin[1], e->origin[2]
	};
	sceGumTranslate(&translation);
/*
	// Scale.
	const ScePspFVector3 scale = {
        e->scale[0],e->scale[1],e->scale[2]
	};
	sceGumScale(&scale);
*/
	// Rotate.
	const ScePspFVector3 rotation = {
		e->angles[ROLL] * (GU_PI / 180.0f),
		-e->angles[PITCH] * (GU_PI / 180.0f),
		e->angles[YAW] * (GU_PI / 180.0f)
	};
	sceGumRotateZYX(&rotation);

	sceGumUpdateMatrix();
}

/*
=============
R_BlendedRotateForEntity

fenix@io.com: model transform interpolation
=============
*/
void R_BlendedRotateForEntity (entity_t *e)
{
    float timepassed;
    float blend;
    vec3_t d;
    int i;

    // positional interpolation

    timepassed = realtime - e->translate_start_time;

    if (e->translate_start_time == 0 || timepassed > 1)
    {
        e->translate_start_time = realtime;
        VectorCopy (e->origin, e->origin1);
        VectorCopy (e->origin, e->origin2);
    }

    if (!VectorCompare (e->origin, e->origin2))
    {
        e->translate_start_time = realtime;
        VectorCopy (e->origin2, e->origin1);
        VectorCopy (e->origin,  e->origin2);
        blend = 0;
    }
    else
    {
        blend =  timepassed / 0.1;
        if (cl.paused || blend > 1)
            blend = 0;
    }

    VectorSubtract (e->origin2, e->origin1, d);

    // Translate.
    const ScePspFVector3 translation = {
    e->origin[0] + (blend * d[0]),
    e->origin[1] + (blend * d[1]),
    e->origin[2] + (blend * d[2])
    };
    sceGumTranslate(&translation);
/*
    // Scale.
    const ScePspFVector3 scale = {
    e->scale[0] + (blend * d[0]),
    e->scale[1] + (blend * d[1]),
    e->scale[2] + (blend * d[2]
    };
    sceGumScale(&scale);
*/
    // orientation interpolation (Euler angles, yuck!)

    timepassed = realtime - e->rotate_start_time;

    if (e->rotate_start_time == 0 || timepassed > 1)
    {
        e->rotate_start_time = realtime;
        VectorCopy (e->angles, e->angles1);
        VectorCopy (e->angles, e->angles2);
    }

    if (!VectorCompare (e->angles, e->angles2))
    {
        e->rotate_start_time = realtime;
        VectorCopy (e->angles2, e->angles1);
        VectorCopy (e->angles,  e->angles2);
        blend = 0;
    }
    else
    {
        blend = timepassed / 0.1;
        if (cl.paused || blend > 1)
            blend = 1;
    }

    VectorSubtract (e->angles2, e->angles1, d);

    // always interpolate along the shortest path
    for (i = 0; i < 3; i++)
    {
        if (d[i] > 180)
        {
            d[i] -= 360;
        }
        else if (d[i] < -180)
        {
            d[i] += 360;
        }
    }

    // Rotate.
    const ScePspFVector3 rotation = {
    (e->angles1[ROLL] + ( blend * d[ROLL])) * (GU_PI / 180.0f),
    (-e->angles1[PITCH] + (-blend * d[PITCH])) * (GU_PI / 180.0f),
    (e->angles1[YAW] + ( blend * d[YAW])) * (GU_PI / 180.0f)
    };
    sceGumRotateZYX(&rotation);

    sceGumUpdateMatrix();
}

/*
=============================================================

  SPRITE MODELS

=============================================================
*/

/*
================
R_GetSpriteFrame
================
*/
mspriteframe_t *R_GetSpriteFrame (entity_t *currententity)
{
	msprite_t		*psprite;
	mspritegroup_t	*pspritegroup;
	mspriteframe_t	*pspriteframe;
	int				i, numframes, frame;
	float			*pintervals, fullinterval, targettime, time;

	psprite = static_cast<msprite_t*>(currententity->model->cache.data);
	frame = currententity->frame;

	if ((frame >= psprite->numframes) || (frame < 0))
	{
		Con_Printf ("R_DrawSprite: no such frame %d\n", frame);
		frame = 0;
	}

	if (psprite->frames[frame].type == SPR_SINGLE)
	{
		pspriteframe = psprite->frames[frame].frameptr;
	}
	else
	{
		pspritegroup = (mspritegroup_t *)psprite->frames[frame].frameptr;
		pintervals = pspritegroup->intervals;
		numframes = pspritegroup->numframes;
		fullinterval = pintervals[numframes-1];

		time = cl.time + currententity->syncbase;

	// when loading in Mod_LoadSpriteGroup, we guaranteed all interval values
	// are positive, so we don't have to worry about division by 0
		targettime = time - ((int)(time / fullinterval)) * fullinterval;

		for (i=0 ; i<(numframes-1) ; i++)
		{
			if (pintervals[i] > targettime)
				break;
		}

		pspriteframe = pspritegroup->frames[i];
	}

	return pspriteframe;
}


/*
=================
R_DrawSpriteModel

=================
*/
void R_DrawSpriteModel (entity_t *e)
{
	vec3_t			point, v_forward, v_right, v_up;
	msprite_t		*psprite;
	mspriteframe_t	*frame;
	float			*s_up, *s_right;
	float			angle, sr, cr;

	// don't even bother culling, because it's just a single polygon without a surface cache
	frame = R_GetSpriteFrame (e);
	psprite = static_cast<msprite_t*>(currententity->model->cache.data);

	switch(psprite->type)
	{
	case SPR_VP_PARALLEL_UPRIGHT: //faces view plane, up is towards the heavens
		v_up[0] = 0;
		v_up[1] = 0;
		v_up[2] = 1;
		s_up = v_up;
		s_right = vright;
		break;
	case SPR_FACING_UPRIGHT: //faces camera origin, up is towards the heavens
		VectorSubtract(currententity->origin, r_origin, v_forward);
		v_forward[2] = 0;
		VectorNormalizeFast(v_forward);
		v_right[0] = v_forward[1];
		v_right[1] = -v_forward[0];
		v_right[2] = 0;
		v_up[0] = 0;
		v_up[1] = 0;
		v_up[2] = 1;
		s_up = v_up;
		s_right = v_right;
		break;
	case SPR_VP_PARALLEL: //faces view plane, up is towards the top of the screen
		s_up = vup;
		s_right = vright;
		break;
	case SPR_ORIENTED: //pitch yaw roll are independent of camera
		AngleVectors (currententity->angles, v_forward, v_right, v_up);
		s_up = v_up;
		s_right = v_right;
		break;
	case SPR_VP_PARALLEL_ORIENTED: //faces view plane, but obeys roll value
		angle = currententity->angles[ROLL] * M_PI_DIV_180;
		sr = sin(angle);
		cr = cos(angle);
		v_right[0] = vright[0] * cr + vup[0] * sr;
		v_up[0] = vright[0] * -sr + vup[0] * cr;
		v_right[1] = vright[1] * cr + vup[1] * sr;
		v_up[1] = vright[1] * -sr + vup[1] * cr;
		v_right[2] = vright[2] * cr + vup[2] * sr;
		v_up[2] = vright[2] * -sr + vup[2] * cr;
		s_up = v_up;
		s_right = v_right;
		break;
	default:
		return;
	}

    additive = false;
	filter = false;

	if (psprite->beamlength == 10) // we use the beam length of sprites, since they are unused by quake anyway.
		additive = true;

	if (psprite->beamlength == 20)
		filter = true;

	// Bind the texture.
	GL_Bind(frame->gl_texturenum);

	sceGuEnable(GU_BLEND);
	sceGuDisable(GU_FOG);

    if (additive)
	{
		sceGuDepthMask(GU_TRUE);
		sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_FIX, 0, 0xFFFFFFFF);
		sceGuTexFunc(GU_TFX_MODULATE , GU_TCC_RGB);
	}
    else if (filter)
	{
		sceGuDepthMask(GU_TRUE);
		sceGuBlendFunc(GU_ADD, GU_FIX, GU_SRC_COLOR, 0, 0);
        sceGuTexFunc(GU_TFX_MODULATE , GU_TCC_RGB);
	}
	else
		sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);

	// Allocate memory for this polygon.
	glvert_t* const	vertices =
		static_cast<glvert_t*>(sceGuGetMemory(sizeof(glvert_t) * 4));

	VectorMA (e->origin, frame->down, s_up, point);
	VectorMA (point, frame->left, s_right, point);

	vertices[0].st[0]	= 0.0f;
	vertices[0].st[1]	= 1.0f;
	vertices[0].xyz[0]	= point[0];
	vertices[0].xyz[1]	= point[1];
	vertices[0].xyz[2]	= point[2];

	VectorMA (e->origin, frame->up, s_up, point);
	VectorMA (point, frame->left, s_right, point);

	vertices[1].st[0]	= 0.0f;
	vertices[1].st[1]	= 0.0f;
	vertices[1].xyz[0]	= point[0];
	vertices[1].xyz[1]	= point[1];
	vertices[1].xyz[2]	= point[2];

	VectorMA (e->origin, frame->up, s_up, point);
	VectorMA (point, frame->right, s_right, point);

	vertices[2].st[0]	= 1.0f;
	vertices[2].st[1]	= 0.0f;
	vertices[2].xyz[0]	= point[0];
	vertices[2].xyz[1]	= point[1];
	vertices[2].xyz[2]	= point[2];

	VectorMA (e->origin, frame->down, s_up, point);
	VectorMA (point, frame->right, s_right, point);

	vertices[3].st[0]	= 1.0f;
	vertices[3].st[1]	= 1.0f;
	vertices[3].xyz[0]	= point[0];
	vertices[3].xyz[1]	= point[1];
	vertices[3].xyz[2]	= point[2];

	// Draw the clipped vertices.
	sceGuDrawArray(
		GU_TRIANGLE_FAN,
		GU_TEXTURE_32BITF | GU_VERTEX_32BITF,
		4, 0, vertices);

	sceGuDepthMask(GU_FALSE);
	sceGuDisable(GU_BLEND);
	sceGuEnable(GU_FOG);
	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);
	sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_ONE_MINUS_SRC_ALPHA, 0, 0);
}

/*
=============================================================

  ALIAS MODELS

=============================================================
*/


#define NUMVERTEXNORMALS	162

extern "C" float	r_avertexnormals[NUMVERTEXNORMALS][3];
float r_avertexnormals[NUMVERTEXNORMALS][3] = {
#include "../anorms.h"
};

vec3_t	shadevector;
//float	shadelight, ambientlight; // LordHavoc: .lit support, removed shadelight and ambientlight

// precalculated dot products for quantized angles
#define SHADEDOT_QUANT 16
float	r_avertexnormal_dots[SHADEDOT_QUANT][256] =
#include "../anorm_dots.h"
;

float	*shadedots = r_avertexnormal_dots[0];

// light lerping - pox@planetquake.com

float *shadedots2 = r_avertexnormal_dots[0];

float lightlerpoffset;

// fenix@io.com: model animation interpolation
int lastposenum0;
int	lastposenum;

// fenix@io.com: model transform interpolation
float old_i_model_transform;

// vertex lighting
float	apitch, ayaw;
vec3_t	vertexlight;

/*
=============
GL_DrawAliasFrame
=============
*/
extern vec3_t lightcolor; // LordHavoc: .lit support
void GL_DrawAliasFrame (aliashdr_t *paliashdr, int posenum, float apitch, float ayaw)
{
	float 	l;
	trivertx_t	*verts;
	int		*order;
	int		count;
	float     r,g,b;
//	vec3_t		l_v;

    lastposenum = posenum;

	verts = (trivertx_t *)((byte *)paliashdr + paliashdr->posedata);
	verts += posenum * paliashdr->poseverts;
	order = (int *)((byte *)paliashdr + paliashdr->commands);

	while (1)
	{
		// get the vertex count and primitive type
		count = *order++;
		if (!count)
			break;		// done
		int prim;
		if (count < 0)
		{
			count = -count;

			if (r_showtris.value)
                prim = GU_LINE_STRIP;
            else
                prim = GU_TRIANGLE_FAN;
		}
		else
		{
		    if (r_showtris.value)
                prim = GU_LINE_STRIP;
            else
                prim = GU_TRIANGLE_STRIP;;
		}

		// Allocate the vertices.
		struct vertex
		{
			float u, v;
			unsigned int color;
			float x, y, z;
		};

		vertex* const out = static_cast<vertex*>(sceGuGetMemory(sizeof(vertex) * count));

		for (int vertex_index = 0; vertex_index < count; ++vertex_index)
		{
			// texture coordinates come from the draw list
			out[vertex_index].u = ((float *)order)[0];
			out[vertex_index].v = ((float *)order)[1];
			order += 2;

			// normals and vertexes come from the frame list

			// LordHavoc: .lit support begin
			//l = shadedots[verts->lightnormalindex] * shadelight; // LordHavoc: original code
			//l = shadedots[verts->lightnormalindex];
			// LordHavoc: .lit support end

			// normals and vertexes come from the frame list
			// blend the light intensity from the two frames together

            if(vlight.value)
                // RIOT - Vertex lighting
                l = VLight_LerpLight(verts->lightnormalindex, verts->lightnormalindex, 1, apitch, ayaw);
            else
            {
                float l1, l2, diff;

                l1 = shadedots[verts->lightnormalindex]; 				// Colored Lighting support
                l2 = shadedots2[verts->lightnormalindex];

                if (l1 != l2)
                {
                    if (l1 > l2)
                    {
                        diff = l1 - l2;
                        diff *= lightlerpoffset;
                        l = l1 - diff;
                    }
                    else
                    {
                        diff = l2 - l1;
                        diff *= lightlerpoffset;
                        l = l1 + diff;
                    }
                }
                else
                    l = l1;
            }

			if (r_model_contrast.value)
				l *= l;

//			if (l > 1.5)
//				l = 1.5;

            r = l * lightcolor[0];
            g = l * lightcolor[1];
            b = l * lightcolor[2];

            // PSP handles colors in bytes, so can't go anything higher then 1.0f)
            // TODO: Use hardware lights to overbright the models, instead of using another pass.
            if (r > 1.0f)
                r = 1.0f;
            if (g > 1.0f)
                g = 1.0f;
            if (b > 1.0f)
                b = 1.0f;

			out[vertex_index].x = verts->v[0];
			out[vertex_index].y = verts->v[1];
			out[vertex_index].z = verts->v[2];

            if (!r_showtris.value)
                out[vertex_index].color = GU_COLOR(r, g, b, currententity->alpha);
            else
                out[vertex_index].color = 0xffffffff;

			++verts;
		}

		if (r_showtris.value)
		{
		    sceGuDisable(GU_TEXTURE_2D);
            sceGumDrawArray(prim, GU_TEXTURE_32BITF | GU_VERTEX_32BITF | GU_COLOR_8888, count, 0, out);
            sceGuEnable(GU_TEXTURE_2D);
		}
        else
            sceGuDrawArray(prim, GU_TEXTURE_32BITF | GU_VERTEX_32BITF | GU_COLOR_8888, count, 0, out);
	}
}

/*
=============
GL_DrawAliasBlendedFrame

fenix@io.com: model animation interpolation
=============
*/
void GL_DrawAliasBlendedFrame (aliashdr_t *paliashdr, int pose1, int pose2, float blend, float apitch, float ayaw)
{
	float       l;
    float       r,g,b;
	trivertx_t* verts1;
	trivertx_t* verts2;
	int*        order;
	int         count, brightness;
	vec3_t      d;

	lastposenum0 = pose1;
	lastposenum  = pose2;

	verts1 = (trivertx_t *)((byte *)paliashdr + paliashdr->posedata);
	verts2 = verts1;

	verts1 += pose1 * paliashdr->poseverts;
	verts2 += pose2 * paliashdr->poseverts;

	order = (int *)((byte *)paliashdr + paliashdr->commands);

	while (1)
	{
       // get the vertex count and primitive type
		int prim;
		count = *order++;

		if (!count) break;

		if (count < 0)
        {
			count = -count;

			if (r_showtris.value)
                prim = GU_LINE_STRIP;
            else
                prim = GU_TRIANGLE_FAN;
		}
		else
		{
		    if (r_showtris.value)
                prim = GU_LINE_STRIP;
            else
                prim = GU_TRIANGLE_STRIP;;
		}

		// Allocate the vertices.
		struct vertex
		{
			float u, v;
			unsigned int color;
			float x, y, z;
		};

		vertex* const out = static_cast<vertex*>(sceGuGetMemory(sizeof(vertex) * count));

		for (int vertex_index = 0; vertex_index < count; ++vertex_index)
		{
			// texture coordinates come from the draw list
			out[vertex_index].u = ((float *)order)[0];
			out[vertex_index].v = ((float *)order)[1];
			order += 2;

			// normals and vertexes come from the frame list
			// blend the light intensity from the two frames together

			if(vlight.value)
			    // RIOT - Vertex lighting
                l = VLight_LerpLight(verts1->lightnormalindex, verts2->lightnormalindex, blend, apitch, ayaw);
            else
            {
                float l1, l2, diff;

                l1 = shadedots[verts1->lightnormalindex]; 				// Colored Lighting support
                l2 = shadedots2[verts1->lightnormalindex];

                if (l1 != l2)
                {
                    if (l1 > l2)
                    {
                        diff = l1 - l2;
                        diff *= lightlerpoffset;
                        l = l1 - diff;
                    }
                    else
                    {
                        diff = l2 - l1;
                        diff *= lightlerpoffset;
                        l = l1 + diff;
                    }
                }
                else
                    l = l1;
            }
			// light contrast - pox@planetquake.com
			if (r_model_contrast.value)
				l *= l;

            r = l * lightcolor[0];
            g = l * lightcolor[1];
            b = l * lightcolor[2];

            if (r > 1.0f)
                r = 1.0f;
            if (g > 1.0f)
                g = 1.0f;
            if (b > 1.0f)
                b = 1.0f;

			VectorSubtract(verts2->v, verts1->v, d);
			// blend the vertex positions from each frame together
            out[vertex_index].x = verts1->v[0] + (blend * d[0]);
            out[vertex_index].y = verts1->v[1] + (blend * d[1]);
            out[vertex_index].z = verts1->v[2] + (blend * d[2]);

            if (!r_showtris.value)
                out[vertex_index].color = GU_COLOR(r, g, b, currententity->alpha);
            else
                out[vertex_index].color = 0xffffffff;

//			byte colorval = ((int) (l*brightness)) & 0xFF;
//			out[vertex_index].color = (colorval << 24) | (colorval << 16) | (colorval << 8) | colorval;

            verts1++;
            verts2++;
		}

        if (r_showtris.value)
		{
		    sceGuDisable(GU_TEXTURE_2D);
            sceGumDrawArray(prim, GU_TEXTURE_32BITF | GU_VERTEX_32BITF | GU_COLOR_8888, count, 0, out);
            sceGuEnable(GU_TEXTURE_2D);
		}
		else
            sceGuDrawArray(prim, GU_TEXTURE_32BITF | GU_VERTEX_32BITF | GU_COLOR_8888, count, 0, out);
	}
}

/*
=============
GL_DrawAliasShadow
=============
*/
extern	vec3_t			lightspot;

void GL_DrawAliasShadow (aliashdr_t *paliashdr, int posenum)
{/*
	float	s, t, l;
	int		i, j;
	int		index;
	trivertx_t	*v, *verts;
	int		list;
	int		*order;
	vec3_t	point;
	float	*normal;
	float	height, lheight;
	int		count;

	lheight = currententity->origin[2] - lightspot[2];

	height = 0;
	verts = (trivertx_t *)((byte *)paliashdr + paliashdr->posedata);
	verts += posenum * paliashdr->poseverts;
	order = (int *)((byte *)paliashdr + paliashdr->commands);

	height = -lheight + 1.0;

	while (1)
	{
		// get the vertex count and primitive type
		count = *order++;
		if (!count)
			break;		// done
		if (count < 0)
		{
			count = -count;
			glBegin (GL_TRIANGLE_FAN);
		}
		else
			glBegin (GL_TRIANGLE_STRIP);

		do
		{
			// texture coordinates come from the draw list
			// (skipped for shadows) glTexCoord2fv ((float *)order);
			order += 2;

			// normals and vertexes come from the frame list
			point[0] = verts->v[0] * paliashdr->scale[0] + paliashdr->scale_origin[0];
			point[1] = verts->v[1] * paliashdr->scale[1] + paliashdr->scale_origin[1];
			point[2] = verts->v[2] * paliashdr->scale[2] + paliashdr->scale_origin[2];

			point[0] -= shadevector[0]*(point[2]+lheight);
			point[1] -= shadevector[1]*(point[2]+lheight);
			point[2] = height;
			glVertex3fv (point);

			verts++;
		} while (--count);

		glEnd ();
	}	*/
}



/*
=================
R_SetupAliasFrame

=================
*/
void R_SetupAliasFrame (int frame, aliashdr_t *paliashdr, float apitch, float ayaw)
{
	int				pose, numposes;
	float			interval;

	if ((frame >= paliashdr->numframes) || (frame < 0))
	{
		Con_DPrintf ("R_AliasSetupFrame: no such frame %d\n", frame);
		frame = 0;
	}

	pose = paliashdr->frames[frame].firstpose;
	numposes = paliashdr->frames[frame].numposes;

	if (numposes > 1)
	{
		interval = paliashdr->frames[frame].interval;
		pose += (int)(cl.time / interval) % numposes;
	}

	GL_DrawAliasFrame (paliashdr, pose, apitch, ayaw);
}

/*
=================
R_SetupAliasBlendedFrame

fenix@io.com: model animation interpolation
=================
*/
void R_SetupAliasBlendedFrame (int frame, aliashdr_t *paliashdr, entity_t* e, float apitch, float ayaw)
{
	int   pose;
	int   numposes;
	float blend;

	if ((frame >= paliashdr->numframes) || (frame < 0))
	{
		Con_DPrintf ("R_AliasSetupFrame: no such frame %d\n", frame);
		frame = 0;
	}

	pose = paliashdr->frames[frame].firstpose;
	numposes = paliashdr->frames[frame].numposes;

	if (numposes > 1)
	{
  		e->frame_interval = paliashdr->frames[frame].interval;
  		pose += (int)(cl.time / e->frame_interval) % numposes;
  	}
    else
    {
		/* One tenth of a second is a good for most Quake animations.
		If the nextthink is longer then the animation is usually meant to pause
		(e.g. check out the shambler magic animation in shambler.qc).  If its
		shorter then things will still be smoothed partly, and the jumps will be
		less noticable because of the shorter time.  So, this is probably a good
		assumption. */
		e->frame_interval = 0.1;
	}

	if (e->pose2 != pose)
	{
		e->frame_start_time = realtime;
		e->pose1 = e->pose2;
		e->pose2 = pose;
		blend = 0;
	}
	else
	{
		blend = (realtime - e->frame_start_time) / e->frame_interval;
	}

	// wierd things start happening if blend passes 1
	if (cl.paused || blend > 1) blend = 1;

	if (blend == 1)
        GL_DrawAliasFrame (paliashdr, pose, apitch, ayaw);
    else
        GL_DrawAliasBlendedFrame (paliashdr, e->pose1, e->pose2, blend, apitch, ayaw);
}

/*
=================
R_DrawAliasModel
=================
*/
void R_DrawAliasModel (entity_t *e)
{
	int			i, j;
	int			lnum;
	vec3_t		dist;
	float		add;
	model_t		*clmodel;
	vec3_t		mins, maxs;
	aliashdr_t	*paliashdr;
	trivertx_t	*verts, *v;
	int			index;
	float		an, s, t;
	int			anim;

//	float		radiusmax = 0.0;

// Tomaz - QC Alpha Begin
/*
    glDisable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
*/
// Tomaz - QC Alpha End

	force_fullbright = false;
	additive = false;
	filter = false;
	alphafunc = false;
    alphafunc2 = false;
	fixlight = false;

	clmodel = currententity->model;

	VectorAdd (currententity->origin, clmodel->mins, mins);
	VectorAdd (currententity->origin, clmodel->maxs, maxs);

	if (R_CullBox (mins, maxs))
		return;

	VectorCopy (currententity->origin, r_entorigin);
	VectorSubtract (r_origin, r_entorigin, modelorg);

	//
	// get lighting information
	//

	// LordHavoc: .lit support begin
	//ambientlight = shadelight = R_LightPoint (currententity->origin); // LordHavoc: original code, removed shadelight and ambientlight
	R_LightPoint(currententity->origin); // LordHavoc: lightcolor is all that matters from this
	// LordHavoc: .lit support end

	// always give the gun some light
	// LordHavoc: .lit support begin
	//if (e == &cl.viewent && ambientlight < 24) // LordHavoc: original code
	//	ambientlight = shadelight = 24; // LordHavoc: original code
/*
	if (e == &cl.viewent)
	{
		if (lightcolor[0] < 24)
			lightcolor[0] = 24;
		if (lightcolor[1] < 24)
			lightcolor[1] = 24;
		if (lightcolor[2] < 24)
			lightcolor[2] = 24;
	}
*/
	// LordHavoc: .lit support end

	for (lnum=0 ; lnum<MAX_DLIGHTS ; lnum++)
	{
		if (cl_dlights[lnum].die >= cl.time)
		{
			VectorSubtract (currententity->origin,
							cl_dlights[lnum].origin,
							dist);
			add = cl_dlights[lnum].radius - Length(dist);

			// LordHavoc: .lit support begin
			/* LordHavoc: original code
			if (add > 0) {
				ambientlight += add;
				//ZOID models should be affected by dlights as well
				shadelight += add;
			}
			*/
			if (add > 0)
			{
				lightcolor[0] += add * cl_dlights[lnum].color[0];
				lightcolor[1] += add * cl_dlights[lnum].color[1];
				lightcolor[2] += add * cl_dlights[lnum].color[2];
			}
			// LordHavoc: .lit support end
		}
	}

	// clamp lighting so it doesn't overbright as much
/*
	if (shadelight > 65)
		shadelight = 65;
	if (ambientlight > 196)
    {
		ambientlight = 196;
		force_fullbright = true;
    }
    else
        force_fullbright = false;
*/
	// ZOID: never allow players to go totally black
//	i = currententity - cl_entities;
//	if (i >= 1 && i<=cl.maxclients /*&& !strcmp (currententity->model->name, "progs/player.mdl") */)
	// LordHavoc: .lit support begin
	//	if (ambientlight < 8) // LordHavoc: original code
	//		ambientlight = shadelight = 8; // LordHavoc: original code

	if (!strcmp (clmodel->name, "progs/eyes.mdl") ||
	    !strcmp (clmodel->name, "progs/flame.mdl") ||
	    !strcmp (clmodel->name, "progs/palmleav.mdl") ||
	    !strcmp (clmodel->name, "progs/k_spike.mdl") ||
	    !strcmp (clmodel->name, "progs/s_spike.mdl") ||
        !strcmp (clmodel->name, "progs/scope.mdl") ||
	    !strcmp (clmodel->name, "progs/spike.mdl"))
	{
	// LordHavoc: .lit support begin
	//	ambientlight = shadelight = 256; // LordHavoc: original code
		lightcolor[0] = lightcolor[1] = lightcolor[2] = 256;
	// LordHavoc: .lit support end
		force_fullbright = true;
	}

	if (!strcmp (clmodel->name, "progs/flame2.mdl") ||
	    !strcmp (clmodel->name, "progs/lavaball.mdl") ||
	    !strcmp (clmodel->name, "progs/bolt.mdl") ||
	    !strcmp (clmodel->name, "progs/bolt2.mdl") ||
	    !strcmp (clmodel->name, "progs/bolt3.mdl") ||
	    !strcmp (clmodel->name, "progs/s_light.mdl") ||
	    !strcmp (clmodel->name, "progs/bullet.mdl") ||
	    !strcmp (clmodel->name, "progs/explode.mdl") ||
	    !strcmp (clmodel->name, "progs/laser.mdl"))
	{
	    lightcolor[0] = lightcolor[1] = lightcolor[2] = 256;
		additive = true;
	}

	if (!strcmp (clmodel->name, "progs/smoke.mdl") ||
	    !strcmp (clmodel->name, "progs/glass1.mdl") ||
	    !strcmp (clmodel->name, "progs/glass2.mdl") ||
	    !strcmp (clmodel->name, "progs/glass3.mdl") ||
	    !strcmp (clmodel->name, "progs/debris.mdl"))
	{
		filter = true;
	}

	if (!strcmp (clmodel->name, "progs/raptor.mdl") ||
        !strcmp (clmodel->name, "progs/5_box_s.mdl") ||
	    !strcmp (clmodel->name, "progs/trex.mdl"))
	{
		alphafunc = true;
	}

	if (!strcmp (clmodel->name, "progs/v_axe.mdl") ||
	    !strcmp (clmodel->name, "progs/v_axea.mdl") ||
	    !strcmp (clmodel->name, "progs/v_bow.mdl") ||
	    !strcmp (clmodel->name, "progs/v_tekbow.mdl") ||
	    !strcmp (clmodel->name, "progs/v_shot.mdl") ||
	    !strcmp (clmodel->name, "progs/v_shot2.mdl") ||
	    !strcmp (clmodel->name, "progs/v_nail.mdl") ||
	    !strcmp (clmodel->name, "progs/v_nail2.mdl") ||
	    !strcmp (clmodel->name, "progs/v_rock.mdl") ||
	    !strcmp (clmodel->name, "progs/v_rock2.mdl") ||
	    !strcmp (clmodel->name, "progs/v_light.mdl") ||
	    !strcmp (clmodel->name, "progs/v_uzi.mdl") ||
	    !strcmp (clmodel->name, "progs/v_sniper.mdl") ||
	    !strcmp (clmodel->name, "progs/crate.mdl") ||
	    !strcmp (clmodel->name, "progs/pc.mdl"))
	{
        if (r_model_brightness.value)
			fixlight = true;
			alphafunc2 = true;
	}

	if (lightcolor[0] < 16)
		lightcolor[0] = 16;
	if (lightcolor[1] < 16)
		lightcolor[1] = 16;
	if (lightcolor[2] < 16)
		lightcolor[2] = 16;

	if (lightcolor[0] > 125 && lightcolor[1] > 125 && lightcolor[2] > 125)
	{
		lightcolor[0] = 125;
		lightcolor[1] = 125;
		lightcolor[2] = 125;
        if (!fixlight)
		    force_fullbright = true;
	}

	// LordHavoc: .lit support end

	// light lerping - pox@planetquake.com
	//shadedots = r_avertexnormal_dots[((int)(e->angles[1] * (SHADEDOT_QUANT / 360.0))) & (SHADEDOT_QUANT - 1)];

    {
		float ang_ceil, ang_floor;

		// add pitch angle so lighting changes when looking up/down (mainly for viewmodel)
    	lightlerpoffset = (e->angles[1]+e->angles[0]) * (SHADEDOT_QUANT / 360.0);

	    ang_ceil = ceil(lightlerpoffset);
    	ang_floor = floor(lightlerpoffset);

    	lightlerpoffset = ang_ceil - lightlerpoffset;

		shadedots = r_avertexnormal_dots[(int)ang_ceil & (SHADEDOT_QUANT - 1)];
		shadedots2 = r_avertexnormal_dots[(int)ang_floor & (SHADEDOT_QUANT - 1)];
    }

	// LordHavoc: .lit support begin
	//shadelight = shadelight / 200.0; // LordHavoc: original code
	VectorScale(lightcolor, 1.0f / 200.0f, lightcolor);
	// LordHavoc: .lit support end

    // light lerping - pox@planetquake.com

//    shadelight = (shadelight + currententity->last_shadelight)/2;
//    currententity->last_shadelight = shadelight;

	an = e->angles[1]/180*M_PI;
	shadevector[0] = cosf(-an);
	shadevector[1] = sinf(-an);
	shadevector[2] = 1;
	VectorNormalize (shadevector);

	//
	// locate the proper data
	//

	paliashdr = (aliashdr_t *)Mod_Extradata (currententity->model);
	c_alias_polys += paliashdr->numtris;

	//
	// draw all the triangles
	//

	sceGumPushMatrix();

	// fenix@io.com: model transform interpolation

	if (r_i_model_transform.value)
		R_BlendedRotateForEntity (e);
	else
		R_RotateForEntity (e);

	const ScePspFVector3 translation = { paliashdr->scale_origin[0], paliashdr->scale_origin[1], paliashdr->scale_origin[2] };
	sceGumTranslate(&translation);

	const ScePspFVector3 scaling = { paliashdr->scale[0], paliashdr->scale[1], paliashdr->scale[2] };
	sceGumScale(&scaling);

	anim = (int)(cl.time*10) & 3;
    GL_Bind(paliashdr->gl_texturenum[currententity->skinnum][anim]);

	// we can't dynamically colormap textures, so they are cached
	// seperately for the players.  Heads are just uncolored.
	if (currententity->colormap != vid.colormap && 0)
	{
		i = currententity - cl_entities;
		if (i >= 1 && i<=cl.maxclients)
		    GL_Bind(playertextures - 1 + i);
	}

	if (force_fullbright)
	{
//	    sceGuDepthMask(GU_TRUE);
	    sceGuEnable(GU_BLEND);
		sceGuEnable(GU_ALPHA_TEST);
		sceGuAlphaFunc(GU_GREATER, 0x88, 0xff);
		sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);
    }
    else if (additive)
	{
	    sceGuEnable(GU_BLEND);
		sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);
		sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_FIX, 0, 0xFFFFFFFF);
	}
    else if (filter)
	{
	    sceGuEnable(GU_BLEND);
		sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);
		sceGuBlendFunc(GU_ADD, GU_FIX, GU_SRC_COLOR, 0, 0);
	}
	else if (alphafunc)
	{
		sceGuEnable(GU_ALPHA_TEST);
		sceGuAlphaFunc(GU_GREATER, 0, 0xff);
		sceGuTexFunc(GU_TFX_MODULATE, GU_TCC_RGBA);
	}
    else if (alphafunc2)
    {
        sceGuEnable(GU_ALPHA_TEST);
        sceGuAlphaFunc(GU_GREATER, 0xaa, 0xff);
        sceGuTexFunc(GU_TFX_MODULATE, GU_TCC_RGBA);
    }
    else
		sceGuTexFunc(GU_TFX_MODULATE, GU_TCC_RGBA);

	sceGuShadeModel(GU_SMOOTH);
	sceGumUpdateMatrix();

    if (r_i_model_animation.value)
        R_SetupAliasBlendedFrame (currententity->frame, paliashdr, currententity, e->angles[0], e->angles[1]);
    else
        R_SetupAliasFrame (currententity->frame, paliashdr, e->angles[0], e->angles[1]);

	if (!force_fullbright && fixlight) //Draws another pass, ouch!
	{
        sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_FIX, 0, 0xFFFFFFFF);
        sceGuEnable(GU_BLEND);

        if (r_i_model_animation.value)
            R_SetupAliasBlendedFrame (currententity->frame, paliashdr, currententity, e->angles[0], e->angles[1]);
        else
            R_SetupAliasFrame (currententity->frame, paliashdr, e->angles[0], e->angles[1]);

        sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_ONE_MINUS_SRC_ALPHA, 0, 0);
        sceGuDisable(GU_BLEND);
	}

	if (force_fullbright)
	{
//	    sceGuDepthMask(GU_FALSE);
		sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_ONE_MINUS_SRC_ALPHA, 0, 0);
		sceGuDisable(GU_BLEND);
	}

	if (additive || filter || alphafunc)
	{
		sceGuBlendFunc(GU_ADD, GU_SRC_ALPHA, GU_ONE_MINUS_SRC_ALPHA, 0, 0);
		sceGuDisable(GU_BLEND);
	}

	sceGuAlphaFunc(GU_GREATER, 0, 0xff);
	sceGuDisable(GU_ALPHA_TEST);
	sceGuTexFunc(GU_TFX_REPLACE , GU_TCC_RGBA);
	sceGuShadeModel(GU_FLAT);

//	sceGuDisable(GU_LIGHTING);

	sceGumPopMatrix();
	sceGumUpdateMatrix();

	if (r_shadows.value)
	{
		sceGumPushMatrix();
		R_RotateForEntity (e);

		/*
		glDisable (GL_TEXTURE_2D);
		glEnable (GL_BLEND);
		glColor4f (0,0,0,0.5);*/
		GL_DrawAliasShadow (paliashdr, lastposenum);
		/*glEnable (GL_TEXTURE_2D);
		glDisable (GL_BLEND);
		glColor4f (1,1,1,1);
		glPopMatrix ();*/

		sceGumPopMatrix();
		sceGumUpdateMatrix();
	}
/*
    // Tomaz - QC Alpha Begin

    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glDisable(GL_BLEND);
    glDisable(GL_ALPHA_TEST);

    // Tomaz - QC Alpha End
*/
}

//==================================================================================

/*
=============
R_DrawEntitiesOnList
=============
*/
void R_DrawEntitiesOnList (void)
{
	int		i;

	if (!r_drawentities.value)
		return;

    // Tomaz - QC Alpha Scale Begin

    if ((!strncmp (currententity->model->name, "progs/bolt", 10)) ||
    (!strncmp (currententity->model->name, "progs/bush", 10)) ||
    (!strncmp (currententity->model->name, "progs/bush1", 11)) ||
    (!strncmp (currententity->model->name, "progs/bush2", 11)) ||
    (!strncmp (currententity->model->name, "progs/bush3", 11)) ||
	(!strncmp (currententity->model->name, "progs/flame", 11)))
    {
        currententity->alpha = 1;
        currententity->scale = 1;
    }

    // Tomaz - QC Alpha Scale End

	// draw sprites seperately, because of alpha blending
	for (i=0 ; i<cl_numvisedicts ; i++)
	{
		currententity = cl_visedicts[i];

		switch (currententity->model->type)
		{
		case mod_alias:
			R_DrawAliasModel (currententity);
			break;
		}

		switch (currententity->model->type)
		{
		case mod_brush:
			R_DrawBrushModel (currententity);
			break;

		default:
			break;
		}
	}

	for (i=0 ; i<cl_numvisedicts ; i++)
	{
		currententity = cl_visedicts[i];

		switch (currententity->model->type)
		{
		case mod_sprite:
			R_DrawSpriteModel (currententity);
			break;
		}
	}
}

// Nehahra - Model_Alpha (Function by FrikaC)
/*
=============
R_DrawTransEntities
=============
*/
/*
void R_DrawTransEntities (void)
{
	// need to draw back to front
	// fixme: this isn't my favorite option
	int		i;
	float bestdist, dist;
	entity_t *bestent;
	vec3_t start, test;

	VectorCopy(r_refdef.vieworg, start);

	if (!r_drawentities.value)
		return;

transgetent:
	bestdist = 0;
	for (i=0 ; i<cl_numvisedicts ; i++)
	{
		currententity = cl_visedicts[i];
		if (currententity->transignore)
			continue;
		if (currententity->transparency == 1 || currententity->transparency ==0)
			continue;

		VectorCopy(currententity->origin, test);
		if (currententity->model->type == mod_brush)
		{
			test[0] += currententity->model->mins[0];
			test[1] += currententity->model->mins[1];
			test[2] += currententity->model->mins[2];
		}
		dist = (((test[0] - start[0]) * (test[0] - start[0])) +
			((test[1] - start[1]) * (test[1] - start[1])) +
			((test[2] - start[2]) * (test[2] - start[2])));

		if (dist > bestdist)
		{
			bestdist = dist;
			bestent = currententity;

		}
	}
	if (bestdist == 0)
		return;
	bestent->transignore = true;

	currententity = bestent;
	switch (currententity->model->type)
	{
	case mod_alias:
        R_DrawAliasModel (currententity);
		break;
	case mod_brush:
		R_DrawBrushModel (currententity);
		break;
	default:
		break;
	}

	goto transgetent;

}
*/
/*
=============
R_DrawViewModel
=============
*/
void R_DrawViewModel (void)
{
	/* LordHavoc: .lit support, this code was completely useless
	float		ambient[4], diffuse[4];
	int			j;
	int			lnum;
	vec3_t		dist;
	float		add;
	dlight_t	*dl;
	int			ambientlight, shadelight;
	*/

	if (!r_drawviewmodel.value || chase_active.value || envmap || !r_drawentities.value || cl.items & IT_INVISIBILITY || cl.stats[STAT_HEALTH] <= 0 || !cl.viewent.model)
		return;

	currententity = &cl.viewent;

    // Tomaz - QC Alpha Scale Begin

    currententity->alpha = cl_entities[cl.viewentity].alpha;
    currententity->scale = cl_entities[cl.viewentity].scale;

    // Tomaz - QC Alpha Scale End

	/* LordHavoc: .lit support, this code was completely useless
	j = R_LightPoint (currententity->origin);

	if (j < 16)
		j = 16;		// allways give some light on gun
	ambientlight = j;
	shadelight = j;

// add dynamic lights
	for (lnum=0 ; lnum<MAX_DLIGHTS ; lnum++)
	{
		dl = &cl_dlights[lnum];
		if (!dl->radius)
			continue;
		if (!dl->radius)
			continue;
		if (dl->die < cl.time)
			continue;

		VectorSubtract (currententity->origin, dl->origin, dist);
		add = dl->radius - Length(dist);
		if (add > 0)
			ambientlight += int(add);
	}

	ambient[0] = ambient[1] = ambient[2] = ambient[3] = (float)ambientlight;// / 128;
	diffuse[0] = diffuse[1] = diffuse[2] = diffuse[3] = (float)shadelight;// / 128;
	*/

	// hack the depth range to prevent view model from poking into walls
	sceGuDepthRange(0, 19660);

         // fenix@io.com: model transform interpolation
         old_i_model_transform = r_i_model_transform.value;
         r_i_model_transform.value = false;
         R_DrawAliasModel (currententity);
         r_i_model_transform.value = old_i_model_transform;

	sceGuDepthRange(0, 65535);

//	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);
//	sceGuDisable(GU_BLEND);
}


/*
============
R_PolyBlend
============
*/
void R_PolyBlend (void)
{
/*	if (!gl_polyblend.value)
		return;*/
	if (!v_blend[3])
		return;
/*
	GL_DisableMultitexture();

	glDisable (GL_ALPHA_TEST);
	glEnable (GL_BLEND);
	glDisable (GL_DEPTH_TEST);
	glDisable (GL_TEXTURE_2D);

    glLoadIdentity ();

    glRotatef (-90,  1, 0, 0);	    // put Z going up
    glRotatef (90,  0, 0, 1);	    // put Z going up

	glColor4fv (v_blend);

	glBegin (GL_QUADS);

	glVertex3f (10, 100, 100);
	glVertex3f (10, -100, 100);
	glVertex3f (10, -100, -100);
	glVertex3f (10, 100, -100);
	glEnd ();

	glDisable (GL_BLEND);
	glEnable (GL_TEXTURE_2D);
	glEnable (GL_ALPHA_TEST);*/
}


static int SignbitsForPlane (mplane_t *out)
{
	int	bits, j;

	// for fast box on planeside test

	bits = 0;
	for (j=0 ; j<3 ; j++)
	{
		if (out->normal[j] < 0)
			bits |= 1<<j;
	}
	return bits;
}


void R_SetFrustum (void)
{
	int		i;

	if (r_refdef.fov_x == 90)
	{
		// front side is visible

		VectorAdd (vpn, vright, frustum[0].normal);
		VectorSubtract (vpn, vright, frustum[1].normal);

		VectorAdd (vpn, vup, frustum[2].normal);
		VectorSubtract (vpn, vup, frustum[3].normal);
	}
	else
	{
		// rotate VPN right by FOV_X/2 degrees
		RotatePointAroundVector( frustum[0].normal, vup, vpn, -(90-r_refdef.fov_x / 2 ) );
		// rotate VPN left by FOV_X/2 degrees
		RotatePointAroundVector( frustum[1].normal, vup, vpn, 90-r_refdef.fov_x / 2 );
		// rotate VPN up by FOV_X/2 degrees
		RotatePointAroundVector( frustum[2].normal, vright, vpn, 90-r_refdef.fov_y / 2 );
		// rotate VPN down by FOV_X/2 degrees
		RotatePointAroundVector( frustum[3].normal, vright, vpn, -( 90 - r_refdef.fov_y / 2 ) );
	}

	for (i=0 ; i<4 ; i++)
	{
		frustum[i].type = PLANE_ANYZ;
		frustum[i].dist = DotProduct (r_origin, frustum[i].normal);
		frustum[i].signbits = SignbitsForPlane (&frustum[i]);
	}
}



/*
===============
R_SetupFrame
===============
*/
void R_SetupFrame (void)
{
	Fog_SetupFrame ();

// don't allow cheats in multiplayer
	if (cl.maxclients > 1)
	{
		Cvar_Set ("r_fullbright", "0");
		Cvar_Set ("r_lightmap", "0");
		Cvar_Set ("r_showtris", "0");
	}

	R_AnimateLight ();

	r_framecount++;

// build the transformation matrix for the given view angles
	VectorCopy (r_refdef.vieworg, r_origin);

	AngleVectors (r_refdef.viewangles, vpn, vright, vup);

// current viewleaf
	r_oldviewleaf = r_viewleaf;
	r_viewleaf = Mod_PointInLeaf (r_origin, cl.worldmodel);

	V_SetContentsColor (r_viewleaf->contents);
	/*
	V_CalcBlend ();
	*/

	r_cache_thrash = qfalse;

	c_brush_polys = 0;
	c_alias_polys = 0;

}


/*
=============
R_SetupGL
=============
*/
void R_SetupGL (void)
{
	float	screenaspect;
	extern	int glwidth, glheight;
	int		x, x2, y2, y, w, h;
	float fovx, fovy; //johnfitz
	int contents; //johnfitz

	//
	// set up viewpoint
	//
	/*
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity ();
	*/
	sceGumMatrixMode(GU_PROJECTION);
	sceGumLoadIdentity();

	x = r_refdef.vrect.x * glwidth/vid.width;
	x2 = (r_refdef.vrect.x + r_refdef.vrect.width) * glwidth/vid.width;
	y = (vid.height-r_refdef.vrect.y) * glheight/vid.height;
	y2 = (vid.height - (r_refdef.vrect.y + r_refdef.vrect.height)) * glheight/vid.height;

	// fudge around because of frac screen scale
	if (x > 0)
		x--;
	if (x2 < glwidth)
		x2++;
	if (y2 < 0)
		y2--;
	if (y < glheight)
		y++;

	w = x2 - x;
	h = y - y2;

	if (envmap)
	{
		x = y2 = 0;
		w = h = 256;
	}

	sceGuViewport(
		glx,
		gly + (glheight / 2) - y2 - (h / 2),
		w,
		h);
	sceGuScissor(x, glheight - y2 - h, x + w, glheight - y2);

    screenaspect = (float)r_refdef.vrect.width/r_refdef.vrect.height;

	//johnfitz -- warp view for underwater
	fovx = screenaspect;
	fovy = r_refdef.fov_y;
//	if (r_waterwarp.value)
//	{
		contents = Mod_PointInLeaf (r_origin, cl.worldmodel)->contents;
		if (contents == CONTENTS_WATER ||
			contents == CONTENTS_SLIME ||
			contents == CONTENTS_LAVA)
		{
			//variance should be a percentage of width, where width = 2 * tan(fov / 2)
			//otherwise the effect is too dramatic at high FOV and too subtle at low FOV
			//what a mess!
			fovx = atan(tan(DEG2RAD(r_refdef.fov_x) / 2) * (0.97 + sin(cl.time * 1) * 0.04)) * 2 / M_PI_DIV_180;
			fovy = atan(tan(DEG2RAD(r_refdef.fov_y) / 2) * (1.03 - sin(cl.time * 2) * 0.04)) * 2 / M_PI_DIV_180;

			//old method where variance was a percentage of fov
			//fovx = r_refdef.fov_x * (0.98 + sin(cl.time * 1.5) * 0.02);
			//fovy = r_refdef.fov_y * (1.02 - sin(cl.time * 1.5) * 0.02);
		}
//	}

	sceGumPerspective(fovy, screenaspect, 4, 4096);

	if (mirror)
	{
		if (mirror_plane->normal[2])
		{
			/*glScalef (1, -1, 1);*/
		}
		else
		{
			/*glScalef (-1, 1, 1);*/
		}
		/*glCullFace(GL_BACK);*/
	}
	else
	{
		/*glCullFace(GL_FRONT);*/
	}
	sceGumUpdateMatrix();

	/*glMatrixMode(GL_MODELVIEW);
    glLoadIdentity ();*/
	sceGumMatrixMode(GU_VIEW);
	sceGumLoadIdentity();

    /*glRotatef (-90,  1, 0, 0);	    // put Z going up*/
	sceGumRotateX(-90 * (GU_PI / 180.0f));

    /*glRotatef (90,  0, 0, 1);	    // put Z going up*/
	sceGumRotateZ(90 * (GU_PI / 180.0f));

    /*glRotatef (-r_refdef.viewangles[2],  1, 0, 0);*/
	sceGumRotateX(-r_refdef.viewangles[2] * (GU_PI / 180.0f));

    /*glRotatef (-r_refdef.viewangles[0],  0, 1, 0);*/
	sceGumRotateY(-r_refdef.viewangles[0] * (GU_PI / 180.0f));

    /*glRotatef (-r_refdef.viewangles[1],  0, 0, 1);*/
	sceGumRotateZ(-r_refdef.viewangles[1] * (GU_PI / 180.0f));

    /*glTranslatef (-r_refdef.vieworg[0],  -r_refdef.vieworg[1],  -r_refdef.vieworg[2]);*/
	const ScePspFVector3 translation = {
		-r_refdef.vieworg[0],
		-r_refdef.vieworg[1],
		-r_refdef.vieworg[2]
	};
	sceGumTranslate(&translation);

	/*glGetFloatv (GL_MODELVIEW_MATRIX, r_world_matrix);*/
	sceGumStoreMatrix(&r_world_matrix);
	sceGumUpdateMatrix();

	sceGumMatrixMode(GU_MODEL);

	clipping::begin_frame();

	//
	// set drawing parms
	//
	/*
	if (gl_cull.value)
	{
		glEnable(GL_CULL_FACE);
	}
	else
	{
		glDisable(GL_CULL_FACE);
	}

	glDisable(GL_BLEND);
	glDisable(GL_ALPHA_TEST);
	glEnable(GL_DEPTH_TEST);*/
	sceGuDisable(GU_BLEND);
	sceGuDisable(GU_ALPHA_TEST);
	sceGuTexFunc(GU_TFX_REPLACE, GU_TCC_RGBA);
}

/*
================
R_RenderScene

r_refdef must be set before the first call
================
*/
void R_RenderScene (void)
{
	R_SetupFrame ();

	R_SetFrustum ();

	R_SetupGL ();

	R_MarkLeaves ();	// done here so we know if we're in water
/*
	// Set up fogging.

    if (r_refdef.fog_end == 0)
		r_refdef.fog_end = -1;
    else
    {
        sceGuEnable ( GU_FOG );
        sceGuFog ( r_refdef.fog_start, r_refdef.fog_end, GU_COLOR( r_refdef.fog_red * 0.01f, r_refdef.fog_green * 0.01f, r_refdef.fog_blue * 0.01f, 1.0f ) );
	}
*/
	Fog_EnableGFog (); //johnfitz

	R_DrawWorld ();		// adds static entities to the list

	S_ExtraUpdate ();	// don't let sound get messed up if going slow

	R_DrawEntitiesOnList ();

	/*GL_DisableMultitexture();*/

	R_RenderDlights ();

	R_DrawWaterSurfaces ();

	R_DrawParticles ();

	R_DrawViewModel ();

	Fog_DisableGFog (); //johnfitz
#ifdef GLTEST
	Test_Draw ();
#endif

}

/*
=============
R_Clear
=============
*/
void R_Clear (void)
{
    if(r_refdef.fog_end > 0 && r_skyfog.value)
    {
        sceGuClear(GU_COLOR_BUFFER_BIT | GU_DEPTH_BUFFER_BIT);
        sceGuClearColor(GU_COLOR(r_refdef.fog_red * 0.01f,r_refdef.fog_green * 0.01f,r_refdef.fog_blue * 0.01f,1.0f));
    }
    else
    {
        sceGuClear(GU_DEPTH_BUFFER_BIT);
//        sceGuClearColor(GU_COLOR(0.25f,0.25f,0.25f,0.5f));
    }

	if (r_mirroralpha.value != 1.0)
	{
		/*
		gldepthmin = 0;
		gldepthmax = 0.5;
		glDepthFunc (GL_LEQUAL);
		*/
	}
	else
	{
		/*
		gldepthmin = 0;
		gldepthmax = 1;
		glDepthFunc (GL_LEQUAL);*/
	}

	/*glDepthRange (gldepthmin, gldepthmax);*/
}

/*
=============
R_Mirror
=============

void R_Mirror (void)
{

	msurface_t	*s;

if(!kurok)
{

	float		d;
	entity_t	*ent;

	if (!mirror)
		return;

	r_base_world_matrix = r_world_matrix;

	d = DotProduct (r_refdef.vieworg, mirror_plane->normal) - mirror_plane->dist;
	VectorMA (r_refdef.vieworg, -2*d, mirror_plane->normal, r_refdef.vieworg);

	d = DotProduct (vpn, mirror_plane->normal);
	VectorMA (vpn, -2*d, mirror_plane->normal, vpn);

	r_refdef.viewangles[0] = -asinf(vpn[2])/M_PI*180;
	r_refdef.viewangles[1] = atan2f(vpn[1], vpn[0])/M_PI*180;
	r_refdef.viewangles[2] = -r_refdef.viewangles[2];

	ent = &cl_entities[cl.viewentity];
	if (cl_numvisedicts < MAX_VISEDICTS)
	{
		cl_visedicts[cl_numvisedicts] = ent;
		cl_numvisedicts++;
	}

	R_RenderScene ();
	R_DrawWaterSurfaces ();

	if (mirror_plane->normal[2])
	{

	}
	else
	{

	}

	s = cl.worldmodel->textures[mirrortexturenum]->texturechain;
	for ( ; s ; s=s->texturechain)
		R_RenderBrushPoly (s);
	cl.worldmodel->textures[mirrortexturenum]->texturechain = NULL;
}
else
{
	s = cl.worldmodel->textures[mirrortexturenum]->texturechain;
	for ( ; s ; s=s->texturechain)
		R_RenderBrushPoly (s);
	cl.worldmodel->textures[mirrortexturenum]->texturechain = NULL;
}

}
*/

/*
======================
R_Fog_f
======================
*/
/*
void R_Fog_f (void)
{
	if (Cmd_Argc () == 1)
	{
		Con_Printf("\"fog\" is \"%f %f %f %f %f\"\n", r_refdef.fog_start, r_refdef.fog_end, r_refdef.fog_red, r_refdef.fog_green, r_refdef.fog_blue);
		return;
	}
	r_refdef.fog_start = atof(Cmd_Argv(1));
	r_refdef.fog_end = atof(Cmd_Argv(2));
	r_refdef.fog_red = atof(Cmd_Argv(3));
	r_refdef.fog_green = atof(Cmd_Argv(4));
	r_refdef.fog_blue = atof(Cmd_Argv(5));
}
*/
/*
================
R_RenderView

r_refdef must be set before the first call
================
*/
void R_RenderView (void)
{
	double	time1, time2;
    int i;

	if (r_norefresh.value)
		return;

	if (!r_worldentity.model || !cl.worldmodel)
		Sys_Error ("R_RenderView: NULL worldmodel");

	time1 = 0;

	if (r_speeds.value)
	{
		/*glFinish ();*/
		time1 = Sys_FloatTime ();
		c_brush_polys = 0;
		c_alias_polys = 0;
	}

	mirror = qfalse;

	R_Clear ();

	// render normal view

	R_RenderScene ();

//	R_DrawTransEntities();

	// render mirror view
//	R_Mirror ();

//	R_PolyBlend ();

	if (r_speeds.value)
	{
		time2 = Sys_FloatTime ();
		Con_Printf ("%3i ms  %4i wpoly %4i epoly\n", (int)((time2-time1)*1000), c_brush_polys, c_alias_polys);
	}
}
