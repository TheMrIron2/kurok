/*
Copyright (C) 1996-1997 Id Software, Inc.

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
#include "quakedef.h"

#ifdef WIN32
#include "winquake.h"
#endif

#ifdef PSP
#include <pspkernel.h>
#include <psputility.h>
#include "net_dgrm.h"

extern cvar_t	accesspoint;
extern cvar_t	r_wateralpha;
extern cvar_t	r_vsync;
extern cvar_t	r_mipmaps;
extern cvar_t	r_mipmaps_bias;
extern cvar_t	in_freelook_analog;
extern cvar_t	in_disable_analog;
extern cvar_t	in_analog_strafe;
extern cvar_t	lookspring;
extern cvar_t	in_zoom_adjust;
extern cvar_t	in_x_axis_adjust;
extern cvar_t	in_y_axis_adjust;
extern cvar_t	r_dithering;
extern cvar_t   r_i_model_animation;
extern cvar_t   t_i_model_transform;
extern cvar_t	show_fps;
extern cvar_t   sv_aim;
extern cvar_t   noexit;


refdef_t	r_refdef;

#endif

extern qboolean bmg_type_changed;

void (*vid_menudrawfn)(void);
void (*vid_menukeyfn)(int key);

enum {m_none, m_main, m_singleplayer, m_load, m_save, m_multiplayer, m_setup, m_net, m_options, m_video, m_keys, m_help, m_quit, m_serialconfig, m_modemconfig, m_lanconfig, m_gameoptions, m_search, m_slist, m_osk} m_state;

void M_Menu_Main_f (void);
	void M_Menu_SinglePlayer_f (void);
		void M_Menu_Load_f (void);
		void M_Menu_Save_f (void);
	void M_Menu_MultiPlayer_f (void);
		void M_Menu_Setup_f (void);
		void M_Menu_Net_f (void);
	void M_Menu_Options_f (void);
		void M_Menu_Keys_f (void);
		void M_Menu_Video_f (void);
	void M_Menu_Help_f (void);
	void M_Menu_Quit_f (void);
void M_Menu_SerialConfig_f (void);
	void M_Menu_ModemConfig_f (void);
void M_Menu_LanConfig_f (void);
void M_Menu_GameOptions_f (void);
void M_Menu_Search_f (void);
void M_Menu_ServerList_f (void);

void M_Main_Draw (void);
	void M_SinglePlayer_Draw (void);
		void M_Load_Draw (void);
		void M_Save_Draw (void);
	void M_MultiPlayer_Draw (void);
		void M_Setup_Draw (void);
		void M_Net_Draw (void);
	void M_Options_Draw (void);
		void M_Keys_Draw (void);
		void M_Video_Draw (void);
	void M_Help_Draw (void);
	void M_Quit_Draw (void);
void M_SerialConfig_Draw (void);
	void M_ModemConfig_Draw (void);
void M_LanConfig_Draw (void);
void M_GameOptions_Draw (void);
void M_Search_Draw (void);
void M_ServerList_Draw (void);

void M_Main_Key (int key);
void M_SinglePlayer_Key (int key);
	void M_Load_Key (int key);
	void M_Save_Key (int key);
void M_MultiPlayer_Key (int key);
	void M_Setup_Key (int key);
	void M_Net_Key (int key);
void M_Options_Key (int key);
	void M_Keys_Key (int key);
	void M_Video_Key (int key);
void M_Help_Key (int key);
void M_Quit_Key (int key);
void M_SerialConfig_Key (int key);
void M_ModemConfig_Key (int key);
void M_LanConfig_Key (int key);
void M_GameOptions_Key (int key);
void M_Search_Key (int key);
void M_ServerList_Key (int key);

void Con_SetOSKActive(qboolean active);
void M_Menu_OSK_f (char *input, char *output, int outlen);


qboolean	m_entersound;		// play after drawing a frame, so caching
								// won't disrupt the sound
qboolean	m_recursiveDraw;

int			m_return_state;

//int         track;

qboolean	m_return_onerror;
char		m_return_reason [32];

#define StartingGame	(m_multiplayer_cursor == 1)
#define JoiningGame		(m_multiplayer_cursor == 0)
#define SerialConfig	(m_net_cursor == 0)
#define DirectConfig	(m_net_cursor == 1)
#define	IPXConfig		(m_net_cursor == 2)
#define	TCPIPConfig		(m_net_cursor == 3)

void M_ConfigureNetSubsystem(void);

/*
================
M_DrawCharacter

Draws one solid graphics character
================
*/
void M_DrawCharacter (int cx, int line, int num)
{
	Draw_Character ( cx + ((vid.width - 320)>>1), line, num);
}

void M_Print (int cx, int cy, char *str)
{
	while (*str)
	{
		M_DrawCharacter (cx, cy, (*str)+128);
		str++;
		cx += 8;
	}
}

void M_PrintWhite (int cx, int cy, char *str)
{
	while (*str)
	{
		M_DrawCharacter (cx, cy, *str);
		str++;
		cx += 8;
	}
}

void M_DrawTransPic (int x, int y, qpic_t *pic)
{
	Draw_TransPic (x + ((vid.width - 320)>>1), y, pic);
}

void M_DrawPic (int x, int y, qpic_t *pic)
{
	Draw_Pic (x + ((vid.width - 320)>>1), y, pic);
}

byte identityTable[256];
byte translationTable[256];

void M_BuildTranslationTable(int top, int bottom)
{
	int		j;
	byte	*dest, *source;

	for (j = 0; j < 256; j++)
		identityTable[j] = j;
	dest = translationTable;
	source = identityTable;
	memcpy (dest, source, 256);

	if (top < 128)	// the artists made some backwards ranges.  sigh.
		memcpy (dest + TOP_RANGE, source + top, 16);
	else
		for (j=0 ; j<16 ; j++)
			dest[TOP_RANGE+j] = source[top+15-j];

	if (bottom < 128)
		memcpy (dest + BOTTOM_RANGE, source + bottom, 16);
	else
		for (j=0 ; j<16 ; j++)
			dest[BOTTOM_RANGE+j] = source[bottom+15-j];
}


void M_DrawTransPicTranslate (int x, int y, qpic_t *pic)
{
	Draw_TransPicTranslate (x + ((vid.width - 320)>>1), y, pic, translationTable);
}


void M_DrawTextBox (int x, int y, int width, int lines)
{
	qpic_t	*p;
	int		cx, cy;
	int		n;

	// draw left side
	cx = x;
	cy = y;
	p = Draw_CachePic ("gfx/box_tl.lmp");
	M_DrawTransPic (cx, cy, p);
	p = Draw_CachePic ("gfx/box_ml.lmp");
	for (n = 0; n < lines; n++)
	{
		cy += 8;
		M_DrawTransPic (cx, cy, p);
	}
	p = Draw_CachePic ("gfx/box_bl.lmp");
	M_DrawTransPic (cx, cy+8, p);

	// draw middle
	cx += 8;
	while (width > 0)
	{
		cy = y;
		p = Draw_CachePic ("gfx/box_tm.lmp");
		M_DrawTransPic (cx, cy, p);
		p = Draw_CachePic ("gfx/box_mm.lmp");
		for (n = 0; n < lines; n++)
		{
			cy += 8;
			if (n == 1)
				p = Draw_CachePic ("gfx/box_mm2.lmp");
			M_DrawTransPic (cx, cy, p);
		}
		p = Draw_CachePic ("gfx/box_bm.lmp");
		M_DrawTransPic (cx, cy+8, p);
		width -= 2;
		cx += 16;
	}

	// draw right side
	cy = y;
	p = Draw_CachePic ("gfx/box_tr.lmp");
	M_DrawTransPic (cx, cy, p);
	p = Draw_CachePic ("gfx/box_mr.lmp");
	for (n = 0; n < lines; n++)
	{
		cy += 8;
		M_DrawTransPic (cx, cy, p);
	}
	p = Draw_CachePic ("gfx/box_br.lmp");
	M_DrawTransPic (cx, cy+8, p);
}

void M_DrawCheckbox (int x, int y, int on)
{
#if 0
	if (on)
		M_DrawCharacter (x, y, 131);
	else
		M_DrawCharacter (x, y, 129);
#endif
	if (on)
		M_Print (x, y, "on");
	else
		M_Print (x, y, "off");
}

//=============================================================================

int m_save_demonum;

/*
================
M_ToggleMenu_f
================
*/
void M_ToggleMenu_f (void)
{
	m_entersound = true;

	if (key_dest == key_menu)
	{
		if (m_state != m_main)
		{
			M_Menu_Main_f ();
			return;
		}
		key_dest = key_game;
		m_state = m_none;
		return;
	}
	if (key_dest == key_console)
	{
		Con_ToggleConsole_f ();
	}
	else
	{
		M_Menu_Main_f ();
	}
}


//=============================================================================
/* MAIN MENU */

int	m_main_cursor;
#define	MAIN_ITEMS	5

void M_Menu_Main_f (void)
{
	if (key_dest != key_menu)
	{
		m_save_demonum = cls.demonum;
		cls.demonum = -1;
	}
	key_dest = key_menu;
	m_state = m_main;
	m_entersound = true;
}


void M_Main_Draw (void)
{
	int		f;
	qpic_t	*p,*b, *s, *m, *o, *h, *q, *t;
	
	if (kurok)
	{
        t = Draw_CachePic ("gfx/menu/title.lmp");
    	M_DrawPic ((320-t->width)/2, 16, t);

        if (m_main_cursor == 0)
            s = Draw_CachePic ("gfx/menu/single_1.lmp");
        else
            s = Draw_CachePic ("gfx/menu/single_0.lmp");
    	M_DrawPic ((320-s->width)/2, 160, s);
    	
        if (m_main_cursor == 1)
            m = Draw_CachePic ("gfx/menu/multi_1.lmp");
        else
            m = Draw_CachePic ("gfx/menu/multi_0.lmp");
    	M_DrawPic ((320-m->width)/2, 176, m);
    	
        if (m_main_cursor == 2)
            o = Draw_CachePic ("gfx/menu/option_1.lmp");
        else
            o = Draw_CachePic ("gfx/menu/option_0.lmp");
    	M_DrawPic ((320-o->width)/2, 192, o);

        if (m_main_cursor == 3)
            h = Draw_CachePic ("gfx/menu/help_1.lmp");
        else
            h = Draw_CachePic ("gfx/menu/help_0.lmp");
    	M_DrawPic ((320-h->width)/2, 208, h);

        if (m_main_cursor == 4)
            q = Draw_CachePic ("gfx/menu/quit_1.lmp");
        else
            q = Draw_CachePic ("gfx/menu/quit_0.lmp");
    	M_DrawPic ((320-q->width)/2, 224, q);

    }
	else
	{
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );

    	p = Draw_CachePic ("gfx/ttl_main.lmp");
    	M_DrawPic ( (320-p->width)/2, 4, p);
    	M_DrawTransPic (72, 32, Draw_CachePic ("gfx/mainmenu.lmp") );

    	f = (int)(host_time * 10)%6;
    	M_DrawTransPic (54, 32 + m_main_cursor * 20,Draw_CachePic( va("gfx/menudot%i.lmp", f+1 ) ) );
    }

        b = Draw_CachePic ("gfx/m_bttns.lmp");
	    M_DrawPic ( (320-b->width)/2, 248, b );
}


void M_Main_Key (int key)
{
	switch (key)
	{
	case K_ESCAPE:
		key_dest = key_game;
		m_state = m_none;
		cls.demonum = m_save_demonum;
		if (cls.demonum != -1 && !cls.demoplayback && cls.state != ca_connected)
			CL_NextDemo ();
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
        if (++m_main_cursor >= MAIN_ITEMS)
            m_main_cursor = 0;
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
        if (--m_main_cursor < 0)
            m_main_cursor = MAIN_ITEMS - 1;
		break;

	case K_ENTER:
		m_entersound = true;

		switch (m_main_cursor)
		{
		case 0:
			M_Menu_SinglePlayer_f ();
			break;

		case 1:
			M_Menu_MultiPlayer_f ();
			break;

		case 2:
			M_Menu_Options_f ();
			break;

		case 3:
			M_Menu_Help_f ();
			break;

//        if(!kurok)
//        {
    		case 4:
	   		  M_Menu_Quit_f ();
	   		  break;
//        }
		}
	}
}

//=============================================================================
/* SINGLE PLAYER MENU */

int	m_singleplayer_cursor;
#define	SINGLEPLAYER_ITEMS	3


void M_Menu_SinglePlayer_f (void)
{
	key_dest = key_menu;
	m_state = m_singleplayer;
	m_entersound = true;
}


void M_SinglePlayer_Draw (void)
{
	int		f;
	qpic_t	*p,*b, *n, *l, *s, *t;

	b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );

	if (kurok)
	{
        t = Draw_CachePic ("gfx/menu/title.lmp");
    	M_DrawPic ((320-t->width)/2, 16, t);

        if (m_singleplayer_cursor == 0)
            n = Draw_CachePic ("gfx/menu/sp/new_1.lmp");
        else
            n = Draw_CachePic ("gfx/menu/sp/new_0.lmp");
    	M_DrawPic ((320-n->width)/2, 160, n);
    	
        if (m_singleplayer_cursor == 1)
            l = Draw_CachePic ("gfx/menu/sp/load_1.lmp");
        else
            l = Draw_CachePic ("gfx/menu/sp/load_0.lmp");
    	M_DrawPic ((320-l->width)/2, 176, l);
    	
        if (m_singleplayer_cursor == 2)
            s = Draw_CachePic ("gfx/menu/sp/save_1.lmp");
        else
            s = Draw_CachePic ("gfx/menu/sp/save_0.lmp");
    	M_DrawPic ((320-s->width)/2, 192, s);
    }
	else
	{
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );
	    p = Draw_CachePic ("gfx/ttl_sgl.lmp");
	    M_DrawPic ( (320-p->width)/2, 4, p);
	    M_DrawTransPic (72, 32, Draw_CachePic ("gfx/sp_menu.lmp") );

	    f = (int)(host_time * 10)%6;

	    M_DrawTransPic (54, 32 + m_singleplayer_cursor * 20,Draw_CachePic( va("gfx/menudot%i.lmp", f+1 ) ) );
	}
}


void M_SinglePlayer_Key (int key)
{
	switch (key)
	{
	case K_ESCAPE:
		M_Menu_Main_f ();
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		if (++m_singleplayer_cursor >= SINGLEPLAYER_ITEMS)
			m_singleplayer_cursor = 0;
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		if (--m_singleplayer_cursor < 0)
			m_singleplayer_cursor = SINGLEPLAYER_ITEMS - 1;
		break;

	case K_ENTER:
		m_entersound = true;

		switch (m_singleplayer_cursor)
		{
		case 0:
			/*
			if (sv.active)
				if (!SCR_ModalMessage("Are you sure you want to\nstart a new game?\n"))
					break;
			*/
			key_dest = key_game;
			if (sv.active)
				Cbuf_AddText ("disconnect\n");
			Cbuf_AddText ("maxplayers 1\n");
			Cbuf_AddText ("map start\n");
			break;

		case 1:
			M_Menu_Load_f ();
			break;

		case 2:
			M_Menu_Save_f ();
			break;
		}
	}
}

//=============================================================================
/* LOAD/SAVE MENU */

int		load_cursor;		// 0 < load_cursor < MAX_SAVEGAMES

#define	MAX_SAVEGAMES		12
char	m_filenames[MAX_SAVEGAMES][SAVEGAME_COMMENT_LENGTH+1];
int		loadable[MAX_SAVEGAMES];

void M_ScanSaves (void)
{
	int		i, j;
	char	name[MAX_OSPATH];
	FILE	*f;
	int		version;

	for (i=0 ; i<MAX_SAVEGAMES ; i++)
	{
		strcpy (m_filenames[i], "--- EMPTY SLOT ---");
		loadable[i] = false;
		sprintf (name, "%s/s%i.sav", com_gamedir, i);
		f = fopen (name, "r");
		if (!f)
			continue;
		fscanf (f, "%i\n", &version);
		fscanf (f, "%79s\n", name);
		strncpy (m_filenames[i], name, sizeof(m_filenames[i])-1);

	// change _ back to space
		for (j=0 ; j<SAVEGAME_COMMENT_LENGTH ; j++)
			if (m_filenames[i][j] == '_')
				m_filenames[i][j] = ' ';
		loadable[i] = true;
		fclose (f);
	}
}

void M_Menu_Load_f (void)
{
	m_entersound = true;
	m_state = m_load;
	key_dest = key_menu;
	M_ScanSaves ();
}


void M_Menu_Save_f (void)
{
	if (!sv.active)
		return;
	if (cl.intermission)
		return;
	if (svs.maxclients != 1)
		return;
	m_entersound = true;
	m_state = m_save;
	key_dest = key_menu;
	M_ScanSaves ();
}


void M_Load_Draw (void)
{
	int		i;
	qpic_t	*p, *b;
 
	b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
 
    if (kurok)
    {
        p = Draw_CachePic ("gfx/menu/sp/load_0.lmp");
        // line cursor
	    M_DrawCharacter (8, 32 + load_cursor*8, 12+((int)(realtime*30)&1));
    }
    else
    {
	    p = Draw_CachePic ("gfx/p_load.lmp");
        // line cursor
	    M_DrawCharacter (8, 32 + load_cursor*8, 12+((int)(realtime*4)&1));
    }
	M_DrawPic ( (320-p->width)/2, 4, p);

	for (i=0 ; i< MAX_SAVEGAMES; i++)
		M_Print (16, 32 + 8*i, m_filenames[i]);


}


void M_Save_Draw (void)
{
	int		i;
	qpic_t	*p, *b;

	b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );

    if (kurok)
    {
        p = Draw_CachePic ("gfx/menu/sp/save_0.lmp");
        // line cursor
	    M_DrawCharacter (8, 32 + load_cursor*8, 12+((int)(realtime*30)&1));
    }
    else
    {
	    p = Draw_CachePic ("gfx/p_save.lmp");
	    // line cursor
	    M_DrawCharacter (8, 32 + load_cursor*8, 12+((int)(realtime*4)&1));
    }
	M_DrawPic ( (320-p->width)/2, 4, p);

	for (i=0 ; i<MAX_SAVEGAMES ; i++)
		M_Print (16, 32 + 8*i, m_filenames[i]);
}


void M_Load_Key (int k)
{
	switch (k)
	{
	case K_ESCAPE:
		M_Menu_SinglePlayer_f ();
		break;

	case K_ENTER:
		S_LocalSound ("misc/menu2.wav");
		if (!loadable[load_cursor])
			return;
		m_state = m_none;
		key_dest = key_game;

	// Host_Loadgame_f can't bring up the loading plaque because too much
	// stack space has been used, so do it now
		SCR_BeginLoadingPlaque ();

	// issue the load command
		Cbuf_AddText (va ("load s%i\n", load_cursor) );
		return;

	case K_UPARROW:
	case K_LEFTARROW:
		S_LocalSound ("misc/menu1.wav");
		load_cursor--;
		if (load_cursor < 0)
			load_cursor = MAX_SAVEGAMES-1;
		break;

	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound ("misc/menu1.wav");
		load_cursor++;
		if (load_cursor >= MAX_SAVEGAMES)
			load_cursor = 0;
		break;
	}
}


void M_Save_Key (int k)
{
	switch (k)
	{
	case K_ESCAPE:
		M_Menu_SinglePlayer_f ();
		break;

	case K_ENTER:
		m_state = m_none;
		key_dest = key_game;
		Cbuf_AddText (va("save s%i\n", load_cursor));
		return;

	case K_UPARROW:
	case K_LEFTARROW:
		S_LocalSound ("misc/menu1.wav");
		load_cursor--;
		if (load_cursor < 0)
			load_cursor = MAX_SAVEGAMES-1;
		break;

	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound ("misc/menu1.wav");
		load_cursor++;
		if (load_cursor >= MAX_SAVEGAMES)
			load_cursor = 0;
		break;
	}
}

//=============================================================================
/* MULTIPLAYER MENU */

int	m_multiplayer_cursor;

#ifdef PSP
#define	MULTIPLAYER_ITEMS	9
#else
#define	MULTIPLAYER_ITEMS	3
#endif

void M_Menu_MultiPlayer_f (void)
{
	key_dest = key_menu;
	m_state = m_multiplayer;
	m_entersound = true;
}


void M_MultiPlayer_Draw (void)
{
	int		f;
	
	qpic_t	*p,*b, *j, *c, *t, *i, *a;
	
    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );

	if (kurok)
	{
//        M_DrawTransPic (72, 32, Draw_CachePic ("gfx/menu/multi_0.lmp") );

        if (m_multiplayer_cursor == 0)
            j = Draw_CachePic ("gfx/menu/mp/join_1.lmp");
        else
            j = Draw_CachePic ("gfx/menu/mp/join_0.lmp");
    	M_DrawPic ((320-j->width)/2, 72, j);
    	
        if (m_multiplayer_cursor == 1)
            c = Draw_CachePic ("gfx/menu/mp/create_1.lmp");
        else
            c = Draw_CachePic ("gfx/menu/mp/create_0.lmp");
    	M_DrawPic ((320-c->width)/2, 88, c);
    	
        if (m_multiplayer_cursor == 2)
            t = Draw_CachePic ("gfx/menu/mp/setup_1.lmp");
        else
            t = Draw_CachePic ("gfx/menu/mp/setup_0.lmp");
    	M_DrawPic ((320-t->width)/2, 104, t);

        if (m_multiplayer_cursor == 3)
            i = Draw_CachePic ("gfx/menu/mp/inf_1.lmp");
        else
            i = Draw_CachePic ("gfx/menu/mp/inf_0.lmp");
    	M_DrawPic ((320-i->width)/2, 128, i);

	    M_DrawCheckbox ((320/2) - ((3*8)/2), 144, tcpipAvailable && !tcpipAdhoc);

        if (m_multiplayer_cursor == 4)
            a = Draw_CachePic ("gfx/menu/mp/adhoc_1.lmp");
        else
            a = Draw_CachePic ("gfx/menu/mp/adhoc_0.lmp");
    	M_DrawPic ((320-a->width)/2, 152, a);
    	
	    M_DrawCheckbox ((320/2) - ((3*8)/2), 168, tcpipAvailable && tcpipAdhoc);
    	
    	if (m_multiplayer_cursor == 5)
	        M_PrintWhite ((320/2) - ((8*8)/2), 184,  "Add Bot");
        else
	        M_Print ((320/2) - ((8*8)/2), 184,  "Add Bot");

    	if (m_multiplayer_cursor == 6)
	        M_PrintWhite ((320/2) - ((12*8)/2), 192,  "Add Team Bot");
        else
	        M_Print ((320/2) - ((12*8)/2), 192,  "Add Team Bot");

	    if (m_multiplayer_cursor == 7)
	        M_PrintWhite ((320/2) - ((10*8)/2), 200,  "Remove Bot");
        else
	        M_Print ((320/2) - ((10*8)/2), 200,  "Remove Bot");

	    if (m_multiplayer_cursor == 8)
	        M_PrintWhite ((320/2) - ((26*8)/2), 216,  "Co-op Player Camera Change");
        else
	        M_Print ((320/2) - ((26*8)/2), 216,  "Co-op Player Camera Change");

	    if (serialAvailable || ipxAvailable || tcpipAvailable)
	    	return;
	    M_PrintWhite ((320/2) - ((27*8)/2), 232, "No Communications Available");
    }
	else
	{
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );
	    
	    p = Draw_CachePic ("gfx/p_multi.lmp");
	    M_DrawPic ( (320-p->width)/2, 4, p);
	    
	    M_DrawTransPic (72, 32, Draw_CachePic ("gfx/mp_menu.lmp") );
	    
	    f = (int)(host_time * 10)%6;

	    M_DrawTransPic (54, 32 + m_multiplayer_cursor * 20,Draw_CachePic( va("gfx/menudot%i.lmp", f+1 ) ) );
	    
	    #ifdef PSP
	    M_Print (72, 97, "Infrastructure ");
	    M_DrawCheckbox (220, 97, tcpipAvailable && !tcpipAdhoc);
	
	    M_Print (72, 117,  "Adhoc          ");
	    M_DrawCheckbox (220, 117, tcpipAvailable && tcpipAdhoc);
        #endif

	    M_Print (72, 137,  "Add Bot        ");
	    M_Print (72, 157,  "Add Team Bot   ");
	    M_Print (72, 177,  "Remove Bot     ");
	    M_Print (72, 197,  "Co-op Player Camera Change");

	    if (serialAvailable || ipxAvailable || tcpipAvailable)
	    	return;
	    M_PrintWhite ((320/2) - ((27*8)/2), 207, "No Communications Available");
	
    }
}


void M_MultiPlayer_Key (int key)
{
	switch (key)
	{
	case K_ESCAPE:
		M_Menu_Main_f ();
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		if (++m_multiplayer_cursor >= MULTIPLAYER_ITEMS)
			m_multiplayer_cursor = 0;
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		if (--m_multiplayer_cursor < 0)
			m_multiplayer_cursor = MULTIPLAYER_ITEMS - 1;
		break;

	case K_ENTER:
		m_entersound = true;
		switch (m_multiplayer_cursor)
		{
		case 0:
			if (serialAvailable || ipxAvailable || tcpipAvailable)
			    M_Menu_LanConfig_f ();
			break;

		case 1:
			if (serialAvailable || ipxAvailable || tcpipAvailable)
			    M_Menu_LanConfig_f ();
			else
			    M_Menu_GameOptions_f ();
			break;

		case 2:
			M_Menu_Setup_f ();
			break;

#ifdef PSP
		case 3:
			Datagram_Shutdown();
			
			tcpipAvailable = !tcpipAvailable;

			if(tcpipAvailable) {
				tcpipAdhoc = false;
				net_driver_to_use = 0;
				Datagram_Init();
			}

			break;

		case 4:
			Datagram_Shutdown();
			
			tcpipAvailable = !tcpipAvailable;

			if(tcpipAvailable) {
				tcpipAdhoc = true;
				net_driver_to_use = 1;
				Datagram_Init();
			}

			break;

	    case 5:	// add bot
		    Cbuf_AddText ("impulse 101\n");
		    break;

	    case 6:	// add team bot
		    Cbuf_AddText ("impulse 100\n");
		    break;

	    case 7:	// remove bot
            Cbuf_AddText ("impulse 102\n");
		    break;
		     
        case 8:	// player camera switch
            if (coop.value)
                Cbuf_AddText ("impulse 103\n");
            break;

#endif
		}
	}
}

//=============================================================================
/* SETUP MENU */

#ifndef PSP
int		setup_cursor = 4;
int		setup_cursor_table[] = {40, 56, 80, 104, 140};
#else
int		setup_cursor = 5;
int		setup_cursor_table[] = {40, 56, 72, 96, 120, 156};
#endif

char	setup_hostname[16];
char	setup_myname[16];
int		setup_oldtop;
int		setup_oldbottom;
int		setup_top;
int		setup_bottom;

#ifndef PSP
#define	NUM_SETUP_CMDS	5
#else
// Define PSP specific variables
#define	NUM_SETUP_CMDS	6
extern int totalAccessPoints;
extern int accessPointNumber[100];
char    setup_accesspoint[64];
#endif

void M_Menu_Setup_f (void)
{
	key_dest = key_menu;
	m_state = m_setup;
	m_entersound = true;
	Q_strcpy(setup_myname, cl_name.string);
	Q_strcpy(setup_hostname, hostname.string);
	setup_top = setup_oldtop = ((int)cl_color.value) >> 4;
	setup_bottom = setup_oldbottom = ((int)cl_color.value) & 15;
	
#ifdef PSP
	if(totalAccessPoints)
	{
	    sceUtilityGetNetParam(accessPointNumber[(int)accesspoint.value], 0, (netData*)setup_accesspoint);
	}
#endif
}

int	M_ColorForMap (int m)
{
	return m < 128 ? m + 8 : m + 8;
}

void M_Setup_Draw (void)
{
	qpic_t	*p,*b;
	int     offset = 0;
	int				top, bottom, tc, bc;

	if (kurok)
	{
	    if (setup_cursor == 0+offset)
		    M_DrawCharacter (168 + 8*strlen(setup_hostname), setup_cursor_table [setup_cursor], 10+((int)(realtime*30)&1));

	    if (setup_cursor == 1+offset)
		    M_DrawCharacter (168 + 8*strlen(setup_myname), setup_cursor_table [setup_cursor], 10+((int)(realtime*30)&1));
    }
	else
	{
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );
	    
	    if (setup_cursor == 0+offset)
		    M_DrawCharacter (168 + 8*strlen(setup_hostname), setup_cursor_table [setup_cursor], 10+((int)(realtime*4)&1));

	    if (setup_cursor == 1+offset)
		    M_DrawCharacter (168 + 8*strlen(setup_myname), setup_cursor_table [setup_cursor], 10+((int)(realtime*4)&1));
    }
    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
	p = Draw_CachePic ("gfx/p_multi.lmp");
	M_DrawPic ( (320-p->width)/2, 4, p);

#ifdef PSP
	offset = 16;

	M_Print (64, 40, "Access Point");
	M_DrawTextBox (160, 32, 16, 1);
	M_Print (168, 40, setup_accesspoint);
#endif

	M_Print (64, 40+offset, "Host name");
	M_DrawTextBox (160, 32+offset, 16, 1);
	M_Print (168, 56, setup_hostname);

	M_Print (64, 56+offset, "Player name");
	M_DrawTextBox (160, 48+offset, 16, 1);
	M_Print (168, 56+offset, setup_myname);

	if(!kurok)
	{
		M_Print (64, 80+offset, "Top color");
		M_Print (64, 104+offset, "Bottom color");
	}
	else
	{
		M_Print (64, 80+offset, "Top color");
		M_Print (64, 104+offset, "Team color");
	}

	M_DrawTextBox (64, 140+offset-8, 14, 1);
	M_Print (72, 140+offset, "Accept Changes");

	p = Draw_CachePic ("gfx/bigbox.lmp");
	M_DrawTransPic (160, 64+offset, p);

		tc = (setup_top & 15)<<4;
		bc = (setup_bottom & 15)<<4;
		top = M_ColorForMap (tc);
		bottom = M_ColorForMap (bc);

	Draw_Fill ( 248, 72+offset, 56, 28, top);
	Draw_Fill ( 248, 72+28+offset, 56, 28, bottom);

	if(!kurok)
	{
		p = Draw_CachePic ("gfx/menuplyr.lmp");
//		M_BuildTranslationTable(setup_top*16, setup_bottom*16);
//		M_DrawTransPicTranslate (172, 72+offset, p);
		M_DrawTransPic (172, 72+offset, p);
	}

	M_DrawCharacter (56, setup_cursor_table [setup_cursor], 12+((int)(realtime*4)&1));

#ifndef PSP
	offset = 0;
#else
	offset = 1;
#endif
}


void M_Setup_Key (int k)
{
	int	l;
	int	offset = 0;

	switch (k)
	{
	case K_ESCAPE:
		M_Menu_MultiPlayer_f ();
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		setup_cursor--;
		if (setup_cursor < 0)
			setup_cursor = NUM_SETUP_CMDS-1;
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		setup_cursor++;
		if (setup_cursor >= NUM_SETUP_CMDS)
			setup_cursor = 0;
		break;

	case K_LEFTARROW:
#ifdef PSP
		if (setup_cursor == 0)
		{
			S_LocalSound ("misc/menu3.wav");
			if(accesspoint.value > 1)
			{
				Cvar_SetValue("accesspoint", accesspoint.value-1);
				sceUtilityGetNetParam(accessPointNumber[(int)accesspoint.value], 0, (netData*)setup_accesspoint);			
			}
		}
		offset = 1;
#endif
		if (setup_cursor < 2+offset)
			return;
		S_LocalSound ("misc/menu3.wav");
		if (setup_cursor == 2+offset)
			setup_top = setup_top - 1;
		if (setup_cursor == 3+offset)
		{
			if(!kurok)
				setup_bottom = setup_bottom - 1;
			else
				setup_bottom = 4;
		}

		break;
	case K_RIGHTARROW:
#ifdef PSP
		if (setup_cursor == 0)
		{
			S_LocalSound ("misc/menu3.wav");
			if(accesspoint.value < totalAccessPoints)
			{
				Cvar_SetValue("accesspoint", accesspoint.value+1);
				sceUtilityGetNetParam(accessPointNumber[(int)accesspoint.value], 0, (netData*)setup_accesspoint);
			}
		}

		offset = 1;
#endif

		if (setup_cursor < 2+offset)
			return;
forward:
		S_LocalSound ("misc/menu3.wav");
		if (setup_cursor == 2+offset)
			setup_top = setup_top + 1;
		if (setup_cursor == 3+offset)
		{
			if(!kurok)
				setup_bottom = setup_bottom + 1;
			else
				setup_bottom = 13;
		}

		break;

	case K_INS:
#ifdef PSP
		offset = 1;
#endif
		if (setup_cursor == 0+offset)
		{
			M_Menu_OSK_f(setup_hostname, setup_hostname, 16);
			break;
		}

		if (setup_cursor == 1+offset)
		{
			M_Menu_OSK_f(setup_myname, setup_myname,16);
			break;
		}
		break;

	case K_ENTER:
#ifdef PSP
		offset = 1;
#endif
		if (setup_cursor == 0+offset || setup_cursor == 1+offset)
			return;

		if (setup_cursor == 2+offset || setup_cursor == 3+offset)
			goto forward;

		if (setup_cursor < 4+offset)
			break;

		// setup_cursor == 4 (OK)
		if (Q_strcmp(cl_name.string, setup_myname) != 0)
			Cbuf_AddText ( va ("name \"%s\"\n", setup_myname) );
		if (Q_strcmp(hostname.string, setup_hostname) != 0)
			Cvar_Set("hostname", setup_hostname);
		if (setup_top != setup_oldtop || setup_bottom != setup_oldbottom)
			Cbuf_AddText( va ("color %i %i\n", setup_top, setup_bottom) );
		m_entersound = true;
		M_Menu_MultiPlayer_f ();
		break;

	case K_BACKSPACE:
#ifdef PSP
		offset = 1;
#endif
		if (setup_cursor == 0+offset)
		{
			if (strlen(setup_hostname))
				setup_hostname[strlen(setup_hostname)-1] = 0;
		}

		if (setup_cursor == 1+offset)
		{
			if (strlen(setup_myname))
				setup_myname[strlen(setup_myname)-1] = 0;
		}
		break;

	default:
		if (k < 32 || k > 127)
			break;

#ifdef PSP
		offset = 1;
#endif

		if (setup_cursor == 0+offset)
		{
			l = strlen(setup_hostname);
			if (l < 15)
			{
				setup_hostname[l+1] = 0;
				setup_hostname[l] = k;
			}
		}
		if (setup_cursor == 1+offset)
		{
			l = strlen(setup_myname);
			if (l < 15)
			{
				setup_myname[l+1] = 0;
				setup_myname[l] = k;
			}
		}
	}

	if (setup_top > 13)
		setup_top = 0;
	if (setup_top < 0)
		setup_top = 13;
	if (setup_bottom > 13)
		setup_bottom = 0;
	if (setup_bottom < 0)
		setup_bottom = 13;
}

//=============================================================================
/* NET MENU */

int	m_net_cursor;
int m_net_items;
int m_net_saveHeight;

char *net_helpMessage [] =
{
/* .........1.........2.... */
  "                        ",
  " Two computers connected",
  "   through two modems.  ",
  "                        ",

  "                        ",
  " Two computers connected",
  " by a null-modem cable. ",
  "                        ",

  " Novell network LANs    ",
  " or Windows 95 DOS-box. ",
  "                        ",
  "(LAN=Local Area Network)",

  " Commonly used to play  ",
  " over the Internet, but ",
  " also used on a Local   ",
  " Area Network.          "
};

void M_Menu_Net_f (void)
{
	key_dest = key_menu;
	m_state = m_net;
	m_entersound = true;
	m_net_items = 4;

	if (m_net_cursor >= m_net_items)
		m_net_cursor = 0;
	m_net_cursor--;
	M_Net_Key (K_DOWNARROW);
}


void M_Net_Draw (void)
{
	int		f;
	qpic_t	*p,*b;
	
	if (!kurok)
	{
		M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );
    }
    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
	p = Draw_CachePic ("gfx/p_multi.lmp");
	M_DrawPic ( (320-p->width)/2, 4, p);

	f = 32;

	if (serialAvailable)
	{
		p = Draw_CachePic ("gfx/netmen1.lmp");
	}
	else
	{
#ifdef WIN32
		p = NULL;
#else
		p = Draw_CachePic ("gfx/dim_modm.lmp");
#endif
	}

	if (p)
		M_DrawTransPic (72, f, p);

	f += 19;

	if (serialAvailable)
	{
		p = Draw_CachePic ("gfx/netmen2.lmp");
	}
	else
	{
#ifdef WIN32
		p = NULL;
#else
		p = Draw_CachePic ("gfx/dim_drct.lmp");
#endif
	}

	if (p)
		M_DrawTransPic (72, f, p);

	f += 19;
	if (ipxAvailable)
		p = Draw_CachePic ("gfx/netmen3.lmp");
	else
		p = Draw_CachePic ("gfx/dim_ipx.lmp");
	M_DrawTransPic (72, f, p);

	f += 19;
	if (tcpipAvailable)
		p = Draw_CachePic ("gfx/netmen4.lmp");
	else
		p = Draw_CachePic ("gfx/dim_tcp.lmp");
	M_DrawTransPic (72, f, p);

	if (m_net_items == 5)	// JDC, could just be removed
	{
		f += 19;
		p = Draw_CachePic ("gfx/netmen5.lmp");
		M_DrawTransPic (72, f, p);
	}

	f = (320-26*8)/2;
	M_DrawTextBox (f, 134, 24, 4);
	f += 8;
	M_Print (f, 142, net_helpMessage[m_net_cursor*4+0]);
	M_Print (f, 150, net_helpMessage[m_net_cursor*4+1]);
	M_Print (f, 158, net_helpMessage[m_net_cursor*4+2]);
	M_Print (f, 166, net_helpMessage[m_net_cursor*4+3]);

	f = (int)(host_time * 10)%6;
	M_DrawTransPic (54, 32 + m_net_cursor * 20,Draw_CachePic( va("gfx/menudot%i.lmp", f+1 ) ) );
}


void M_Net_Key (int k)
{
again:
	switch (k)
	{
	case K_ESCAPE:
		M_Menu_MultiPlayer_f ();
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		if (++m_net_cursor >= m_net_items)
			m_net_cursor = 0;
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		if (--m_net_cursor < 0)
			m_net_cursor = m_net_items - 1;
		break;

	case K_ENTER:
		m_entersound = true;

		switch (m_net_cursor)
		{
		case 0:
			M_Menu_SerialConfig_f ();
			break;

		case 1:
			M_Menu_SerialConfig_f ();
			break;

		case 2:
			M_Menu_LanConfig_f ();
			break;

		case 3:
			M_Menu_LanConfig_f ();
			break;

		case 4:
// multiprotocol
			break;
		}
	}

	if (m_net_cursor == 0 && !serialAvailable)
		goto again;
	if (m_net_cursor == 1 && !serialAvailable)
		goto again;
	if (m_net_cursor == 2 && !ipxAvailable)
		goto again;
	if (m_net_cursor == 3 && !tcpipAvailable)
		goto again;
}

//=============================================================================
/* OPTIONS MENU */
#define	SLIDER_RANGE	10
#define NUM_SUBMENU 4
#define KNUM_SUBMENU 4
#ifdef PSP_HARDWARE_VIDEO
enum 
{
	OPT_CUSTOMIZE = 0,
	OPT_CONSOLE,
	OPT_GAP,
	OPT_DEFAULTS1,
//	OPT_DEFAULTS2,
//	OPT_DEFAULTS3,
	OPT_GAP1,
	OPT_SUBMENU,
	OPTIONS_ITEMS
};
enum 
{
	OPT_SUBMENU_0 = OPT_SUBMENU,
    OPT_GAP_0,
	OPT_MOUSELOOK,
	OPT_MOUSESTRAFE,
	OPT_AUTOCENTER,
    OPT_GAP_0_1,
	OPT_IN_SPEED,
	OPT_IN_X_ADJUST,
	OPT_IN_Y_ADJUST,
	OPT_IN_TOLERANCE,
	OPT_IN_ACCELERATION,
    OPT_GAP_0_2,
	OPT_ALWAYRUN,
	OPT_INVMOUSE,
	OPT_NOMOUSE,
    OPTIONS_ITEMS_0
};
enum 
{
	OPT_SUBMENU_1 = OPT_SUBMENU,
    OPT_GAP_1,
	OPT_GAMMA,
	OPT_MAXFPS,
    OPT_GAP_1_1,
    OPT_DYNAMIC,
	OPT_MODEL_BRIGHTNESS,
	OPT_SIMPLE_PART,
	OPT_MIPMAPS,
	OPT_ANTIALIAS,
	OPT_VSYNC,
    OPT_FPS,
    OPTIONS_ITEMS_1
};

enum 
{
	OPT_SUBMENU_2 = OPT_SUBMENU,
    OPT_GAP_2,
//    OPT_MUSICTRACK,
	OPT_MUSICVOL,
	OPT_SNDVOL,
    OPT_GAP_2_1,
	OPT_MUSICTYPE,
    OPTIONS_ITEMS_2
};

enum 
{
	OPT_SUBMENU_3 = OPT_SUBMENU,
    OPT_GAP_3,
	OPT_FOG_START,
	OPT_FOG_END,
	OPT_FOG_RED,
	OPT_FOG_GREEN,
	OPT_FOG_BLUE,
    OPT_GAP_3_1,
    OPT_SCRSIZE,
	OPT_WATERTRANS,
	OPT_MIPMAP_BIAS,
	OPT_TEX_SCALEDOWN,
	OPT_SMOOTH_ANIMS,
	OPT_SMOOTH_MOVEMENT,
    OPTIONS_ITEMS_3
};

#else
enum 
{
	OPT_CUSTOMIZE = 0,
	OPT_CONSOLE,
	OPT_DEFAULTS,
	OPT_SUBMENU,
	OPTIONS_ITEMS
};
enum 
{
	OPT_SUBMENU_0 = OPT_SUBMENU,
    //OPT_GAP_0,
	OPT_SCRSIZE,	
	OPT_GAMMA,		
	OPT_VSYNC,	
    //OPT_GAP_0_1,
	OPT_MUSICTYPE,
	OPT_MUSICVOL,
	OPT_SNDVOL,
    OPTIONS_ITEMS_0
};
enum 
{
	OPT_SUBMENU_1 = OPT_SUBMENU,
    //OPT_GAP_1,
	OPT_ALWAYRUN,
	OPT_IN_SPEED,
	OPT_IN_TOLERANCE,
	OPT_IN_ACCELERATION,	
	OPT_INVMOUSE,	
	OPT_NOMOUSE,
	OPT_MOUSELOOK,
	OPT_MOUSESTAFE,
	OPT_IN_X_ADJUST,
	OPT_IN_Y_ADJUST,	
    OPTIONS_ITEMS_1
};
#endif

int	options_cursor;
int m_submenu = 0;

void M_Menu_Options_f (void)
{
	key_dest = key_menu;
	m_state = m_options;
	m_entersound = true;
}

#ifdef PSP_HARDWARE_VIDEO
extern int changeMp3Volume;
#endif

void M_AdjustSliders (int dir)
{
	S_LocalSound ("misc/menu3.wav");

	switch (options_cursor)
	{
		case OPT_SUBMENU:
	        m_submenu += dir;
	        if (kurok)
	        {
	            if (m_submenu > KNUM_SUBMENU-1)
	        	    m_submenu = 0;
	            else if (m_submenu < 0)
	        	    m_submenu = KNUM_SUBMENU-1;
	            break;
            }
	        else
	        {
	            if (m_submenu > NUM_SUBMENU-1)
	        	    m_submenu = 0;
	            else if (m_submenu < 0)
	        	    m_submenu = NUM_SUBMENU-1;
	            break;
            }
	}

    if (m_submenu == 0)
    {
    	switch (options_cursor)
        {
       		case OPT_IN_SPEED:	// mouse speed
				in_sensitivity.value += dir * 1;
				if (in_sensitivity.value < 1)
					in_sensitivity.value = 1;
				if (in_sensitivity.value > 33)
					in_sensitivity.value = 33;
				Cvar_SetValue ("sensitivity", in_sensitivity.value);
				break;

       		case OPT_IN_TOLERANCE:	// mouse tolerance
				in_tolerance.value += dir * 0.01;
				if (in_tolerance.value < 0)
					in_tolerance.value = 0;
				if (in_tolerance.value > 0.5)
					in_tolerance.value = 0.5;
				Cvar_SetValue ("tolerance", in_tolerance.value);
				break;

       		case OPT_IN_ACCELERATION:	// mouse tolerance
				in_acceleration.value -= dir * 0.25;
				if (in_acceleration.value < 0.5)
					in_acceleration.value = 0.5;
				if (in_acceleration.value > 2)
					in_acceleration.value = 2;
				Cvar_SetValue ("acceleration", in_acceleration.value);
				break;
				
			case OPT_IN_X_ADJUST:	
				in_x_axis_adjust.value += dir * 1;
				if (in_x_axis_adjust.value < 1)
					in_x_axis_adjust.value = 1;
				if (in_x_axis_adjust.value > 11)
					in_x_axis_adjust.value = 11;
				Cvar_SetValue ("in_x_axis_adjust", in_x_axis_adjust.value);
				break;
				
			case OPT_IN_Y_ADJUST:	
				in_y_axis_adjust.value += dir * 1;
				if (in_y_axis_adjust.value < 1)
					in_y_axis_adjust.value = 1;
				if (in_y_axis_adjust.value > 11)
					in_y_axis_adjust.value = 11;
				Cvar_SetValue ("in_y_axis_adjust", in_y_axis_adjust.value);
				break;
/*
			case OPT_IN_ZOOM_ADJUST; // zoom speed
				in_zoom_adjust.value += dir * 1;
				if (in_zoom_adjust.value < 1)
					in_zoom_adjust.value = 1;
				if (in_zoom_adjust.value > 33)
					in_zoom_adjust.value = 33;
				Cvar_SetValue ("in_zoom_adjust", in_zoom_adjust.value);
				break;
*/
			case OPT_INVMOUSE:	// invert mouse
				Cvar_SetValue ("m_pitch", -m_pitch.value);
				break;
				
			case OPT_NOMOUSE:	// disable mouse
				Cvar_SetValue ("in_disable_analog", !in_disable_analog.value);
				break;

			case OPT_AUTOCENTER: // auto center looking for digital keys
				Cvar_SetValue ("lookcenter", !lookcenter.value);
				break;

			case OPT_MOUSESTRAFE:
				Cvar_SetValue ("in_analog_strafe", !in_analog_strafe.value);
				break;
				
			case OPT_MOUSELOOK:
				Cvar_SetValue ("in_freelook_analog", !in_freelook_analog.value);
				break;

			case OPT_ALWAYRUN:	// allways run
	            if (kurok)
	            {
		            if (cl_forwardspeed.value > 150)
	                {
        			    Cvar_SetValue ("cl_forwardspeed", 150);
           			    Cvar_SetValue ("cl_backspeed", 150);
        	    	}
        		    else
        		    {
        		    	Cvar_SetValue ("cl_forwardspeed", 200);
        		    	Cvar_SetValue ("cl_backspeed", 200);
        	    	}
                }
            	else
            	{
            		if (cl_forwardspeed.value > 200)
            		{
            			Cvar_SetValue ("cl_forwardspeed", 200);
            			Cvar_SetValue ("cl_backspeed", 200);
            		}
            		else
            		{
        	    		Cvar_SetValue ("cl_forwardspeed", 400);
        	    		Cvar_SetValue ("cl_backspeed", 400);
                    }
        		}
        		break;
        }
    }
    else if (m_submenu == 1)
    {
       	switch (options_cursor)
        {				
			case OPT_GAMMA:	// gamma
				v_gamma.value -= dir * 0.05;
				if (v_gamma.value < 0.5)
					v_gamma.value = 0.5;
				if (v_gamma.value > 1)
					v_gamma.value = 1;
				Cvar_SetValue ("gamma", v_gamma.value);
				break;
				
	        case OPT_MAXFPS:	// frame rate controller
		        max_fps.value += dir * 5;
		        if (max_fps.value < 30)
			        max_fps.value = 30;
		        if (max_fps.value > 65)
			        max_fps.value = 65;
                Cvar_SetValue ("max_fps", max_fps.value);
		        break;
		        
			case OPT_VSYNC:	
				Cvar_SetValue ("r_vsync", !r_vsync.value);
				break;

			case OPT_MODEL_BRIGHTNESS:
				Cvar_SetValue ("r_model_brightness", !r_model_brightness.value);
				break;

			case OPT_FPS:
				Cvar_SetValue ("show_fps", !show_fps.value);
				break;

            case OPT_DYNAMIC:
				Cvar_SetValue ("r_dynamic", !r_dynamic.value);
				break;
#ifdef PSP_HARDWARE_VIDEO
			case OPT_SIMPLE_PART:	
				Cvar_SetValue ("r_particles_simple", !r_particles_simple.value);
				break;
				
			case OPT_ANTIALIAS:	
				Cvar_SetValue ("r_antialias", !r_antialias.value);
				break;

			case OPT_MIPMAPS:
				Cvar_SetValue ("r_mipmaps", !r_mipmaps.value);
				break;
#endif	
        }	
    }
    else if (m_submenu == 2)
    {
       	switch (options_cursor)
        {
/*
			case OPT_MUSICTRACK:	
				track += dir * 1;
				if (track < 1)
					track = 1;
				if (track > 13)
					track = 13;
				Cvar_SetValue ("cd play", track);
				break;
*/
			case OPT_MUSICVOL:	// music volume
#ifdef WIN32
				bgmvolume.value += dir * 1.0;
#else
				bgmvolume.value += dir * 0.1;
#endif
				if (bgmvolume.value < 0)
					bgmvolume.value = 0;
				if (bgmvolume.value > 1)
					bgmvolume.value = 1;
				Cvar_SetValue ("bgmvolume", bgmvolume.value);
#ifdef PSP_MP3HARDWARE_MP3LIB
		        changeMp3Volume = 1;
#endif
				break;
				
			case OPT_SNDVOL:	// sfx volume
				volume.value += dir * 0.1;
				if (volume.value < 0)
					volume.value = 0;
				if (volume.value > 1)
					volume.value = 1;
				Cvar_SetValue ("volume", volume.value);
				break;
				
			case OPT_MUSICTYPE: // bgm type
				if (strcmpi(bgmtype.string,"cd") == 0)
				{
						Cvar_Set("bgmtype","none");
						bmg_type_changed = true;
				}
				else
				{
						Cvar_Set("bgmtype","cd");
						bmg_type_changed = true;
				}
				break;
        }	
    }
    else if (m_submenu == 3)
    {
       	switch (options_cursor)
        {
			case OPT_FOG_START:	// Fog start distance from viewpoint
				r_refdef.fog_start += dir * 100;
				if (r_refdef.fog_start < -5000)
					r_refdef.fog_start = -5000;
				if (r_refdef.fog_start > 5000)
					r_refdef.fog_start = 5000;
				break;

			case OPT_FOG_END:	// Fog end distance from viewpoint
				r_refdef.fog_end += dir * 100;
				if (r_refdef.fog_end < -5000)
					r_refdef.fog_end = -5000;
				if (r_refdef.fog_end > 5000)
					r_refdef.fog_end = 5000;
				break;

			case OPT_FOG_RED:	// Fog red
				r_refdef.fog_red += dir * 5;
				if (r_refdef.fog_red < 0)
					r_refdef.fog_red = 0;
				if (r_refdef.fog_red > 100)
					r_refdef.fog_red = 100;
				break;

			case OPT_FOG_GREEN:	// Fog green
				r_refdef.fog_green += dir * 5;
				if (r_refdef.fog_green < 0)
					r_refdef.fog_green = 0;
				if (r_refdef.fog_green > 100)
					r_refdef.fog_green = 100;
				break;

			case OPT_FOG_BLUE:	// Fog blue
				r_refdef.fog_blue += dir * 5;
				if (r_refdef.fog_blue < 0)
					r_refdef.fog_blue = 0;
				if (r_refdef.fog_blue > 100)
					r_refdef.fog_blue = 100;
				break;

			case OPT_SCRSIZE:	// screen size
				scr_viewsize.value += dir * 10;
				if (scr_viewsize.value < 30)
					scr_viewsize.value = 30;
				if (scr_viewsize.value > 130)
					scr_viewsize.value = 130;
				Cvar_SetValue ("viewsize", scr_viewsize.value);
				break;

			case OPT_WATERTRANS:	// wateralpha
				r_wateralpha.value += dir * 0.1;
				if (r_wateralpha.value < 0)
					r_wateralpha.value = 0;
				if (r_wateralpha.value > 1)
					r_wateralpha.value = 1;
				
				Cvar_SetValue ("r_wateralpha", r_wateralpha.value);
				break;

			case OPT_MIPMAP_BIAS:	// mipmapping bais
				r_mipmaps_bias.value += dir * 0.5;
				if (r_mipmaps_bias.value < -10)
					r_mipmaps_bias.value = -10;
				if (r_mipmaps_bias.value > 0)
					r_mipmaps_bias.value = 0;

				Cvar_SetValue ("r_mipmaps_bias", r_mipmaps_bias.value);
				break;

			case OPT_TEX_SCALEDOWN:	
				Cvar_SetValue ("r_tex_scale_down", !r_tex_scale_down.value);
				break;

			case OPT_SMOOTH_ANIMS:
				Cvar_SetValue ("r_i_model_animation", !r_i_model_animation.value);
				break;

			case OPT_SMOOTH_MOVEMENT:
				Cvar_SetValue ("r_i_model_transform", !r_i_model_transform.value);
				break;
        }	
    }
}


void M_DrawSlider (int x, int y, float range)
{
	int	i;

	if (range < 0)
		range = 0;
	if (range > 1)
		range = 1;
	M_DrawCharacter (x-8, y, 128);
	for (i=0 ; i<SLIDER_RANGE ; i++)
		M_DrawCharacter (x + i*8, y, 129);
	M_DrawCharacter (x+i*8, y, 130);
	M_DrawCharacter (x + (SLIDER_RANGE-1)*8 * range, y, 131);
}

void M_Options_Draw (void)
{
	float		r;
	qpic_t	*p,*b;
	int offset = 32;
	
	if (kurok)
	{
        offset = 64;
	    p = Draw_CachePic ("gfx/menu/option_0.lmp");
	    M_DrawPic ( (320-p->width)/2, 36, p);
	    
	    // Cursor
	    M_DrawCharacter (200, offset + options_cursor*8, 12+((int)(realtime*30)&1));
    }
	else
	{
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );
	    p = Draw_CachePic ("gfx/p_option.lmp");
	    M_DrawPic ( (320-p->width)/2, 4, p);
	    
	    // Cursor
	    M_DrawCharacter (200, offset + options_cursor*8, 12+((int)(realtime*4)&1));
    }

    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );

	M_Print (16, offset+(OPT_CUSTOMIZE*8), "     Customize Buttons");
	M_Print (16, offset+(OPT_CONSOLE*8),   "         Go to console");
	
	M_Print (16, offset+(OPT_DEFAULTS1*8), "         Load defaults");
//	M_Print (16, offset+(OPT_DEFAULTS2*8), "     'Golden' defaults");
//	M_Print (16, offset+(OPT_DEFAULTS3*8), "    'Digital' defaults");
	
	switch (m_submenu)
    {
        case 0:

		#ifdef PSP
			M_Print (16, offset+(OPT_IN_SPEED*8), 		 "          Analog Speed");
		#else
			M_Print (16, offset+(OPT_IN_SPEED*8), 		 "           Mouse Speed");
		#endif
			r = (in_sensitivity.value - 1)/33;
			M_DrawSlider (220, offset+(OPT_IN_SPEED*8), r);

			M_Print (16, offset+(OPT_IN_ACCELERATION*8), "   Analog Acceleration");
			r = 1.0f -((in_acceleration.value - 0.5f)/1.5f);
			M_DrawSlider (220, offset+(OPT_IN_ACCELERATION*8), r);

			M_Print (16, offset+(OPT_IN_TOLERANCE*8), 	 "     Analog Tolerance");
			r = (in_tolerance.value )/1.0f;
			M_DrawSlider (220, offset+(OPT_IN_TOLERANCE*8), r);
		
			M_Print (16, offset+(OPT_IN_X_ADJUST*8), 	 "   Analog Speed X Axis");
			r = (in_x_axis_adjust.value - 1)/10;
			M_DrawSlider (220, offset+(OPT_IN_X_ADJUST*8), r);

			M_Print (16, offset+(OPT_IN_Y_ADJUST*8), 	 "   Analog Speed Y Axis");
			r = (in_y_axis_adjust.value - 1)/10;
			M_DrawSlider (220, offset+(OPT_IN_Y_ADJUST*8), r);
/*
			M_Print (16, offset+(OPT_IN_ZOOM_ADJUST*8),  "            Zoom Speed");
			r = (in_zoom_adjust.value - 1)/33;
			M_DrawSlider (220, offset+(OPT_IN_ZOOM_ADJUST*8), r);
*/
		#ifdef PSP
			M_Print (16, offset+(OPT_INVMOUSE*8),        "         Invert Analog");
		#else
			M_Print (16, offset+(OPT_INVMOUSE*8),        "          Invert Mouse");
		#endif
			M_DrawCheckbox (220, offset+(OPT_INVMOUSE*8), m_pitch.value < 0);
		
			M_Print (16, offset+(OPT_MOUSELOOK*8),       "           Analog Mode");
			if (in_freelook_analog.value == 1)
				M_Print (220, offset+(OPT_MOUSELOOK*8), "Look");
			else
				M_Print (220, offset+(OPT_MOUSELOOK*8), "Move");
//			M_DrawCheckbox (220, offset+(OPT_MOUSELOOK*8), in_freelook_analog.value);

			M_Print (16, offset+(OPT_NOMOUSE*8),         "        Disable Analog");
			M_DrawCheckbox (220, offset+(OPT_NOMOUSE*8), in_disable_analog.value );

			M_Print (16, offset+(OPT_AUTOCENTER*8),      "Autocenter Button Look");
			M_DrawCheckbox (220, offset+(OPT_AUTOCENTER*8), !lookcenter.value );

			M_Print (16, offset+(OPT_MOUSESTRAFE*8),		 "       Analog Strafing");
			M_DrawCheckbox (220, offset+(OPT_MOUSESTRAFE*8), in_analog_strafe.value );

        	M_Print (16, offset+(OPT_ALWAYRUN*8),        "            Always Run");
        	if (kurok)
	            M_DrawCheckbox (220, offset+(OPT_ALWAYRUN*8), cl_forwardspeed.value > 150);
            else
	            M_DrawCheckbox (220, offset+(OPT_ALWAYRUN*8), cl_forwardspeed.value > 200);

			break;

        case 1:
		
			M_Print (16, offset+(OPT_GAMMA*8), 	   "            Brightness");
			r = (1.0 - v_gamma.value) / 0.5;
			M_DrawSlider (220, offset+(OPT_GAMMA*8), r);

			M_Print (16, offset+(OPT_VSYNC*8),     "                 VSync");
			M_DrawCheckbox (220, offset+(OPT_VSYNC*8), r_vsync.value);

			M_Print (16, offset+(OPT_DYNAMIC*8),   "      Dynamic Lighting");
			M_DrawCheckbox (220, offset+(OPT_DYNAMIC*8), r_dynamic.value);

			M_Print (16, offset+(OPT_MODEL_BRIGHTNESS*8),	"       Brighter Models");
			M_DrawCheckbox (220, offset+(OPT_MODEL_BRIGHTNESS*8), r_model_brightness.value);

	        M_Print (16, offset+(OPT_MAXFPS*8),    "    Maximum Frame Rate");
	        r = (max_fps.value - 30) / (65 - 30);
	        M_DrawSlider (220, offset+(OPT_MAXFPS*8), r);

			M_Print (16, offset+(OPT_FPS*8), "    Display Frame rate");
			M_DrawCheckbox (220, offset+(OPT_FPS*8), show_fps.value);

		#ifdef PSP_HARDWARE_VIDEO
		
			M_Print (16, offset+(OPT_SIMPLE_PART*8),    "      Simple Particles");
			M_DrawCheckbox (220, offset+(OPT_SIMPLE_PART*8), r_particles_simple.value);

			M_Print (16, offset+(OPT_ANTIALIAS*8),      "         Anti-Aliasing");
			M_DrawCheckbox (220, offset+(OPT_ANTIALIAS*8), r_antialias.value);

			M_Print (16, offset+(OPT_MIPMAPS*8),        "            MipMapping");
			M_DrawCheckbox (220, offset+(OPT_MIPMAPS*8), r_mipmaps.value);
		#endif

            break;

        case 2:
/*
			M_Print (16, offset+(OPT_MUSICTRACK*8), 	 "           Music Track");
			r = (track - 1)/13;
			M_DrawSlider (220, offset+(OPT_MUSICTRACK*8), r);
*/
			M_Print (16, offset+(OPT_MUSICVOL*8), "          Music Volume");
			r = bgmvolume.value;
			M_DrawSlider (220, offset+(OPT_MUSICVOL*8), r);
		
			M_Print (16, offset+(OPT_SNDVOL*8),   "          Sound Volume");
			r = volume.value;
			M_DrawSlider (220, offset+(OPT_SNDVOL*8), r);

			M_Print (16, offset+(OPT_MUSICTYPE*8),"          MP3 Playback");
			if (strcmpi(bgmtype.string,"cd") == 0)
				M_Print (220, offset+(OPT_MUSICTYPE*8), "On");
			else
				M_Print (220, offset+(OPT_MUSICTYPE*8), "Off");

            break;

        case 3:

			M_Print (16, offset+(OPT_FOG_START*8),			"    Fog start distance");
			r = (r_refdef.fog_start - (-5000)) / ((5000) - (-5000));
			M_DrawSlider (220, offset+(OPT_FOG_START*8), r);
		
			M_Print (16, offset+(OPT_FOG_END*8),			"      Fog end distance");
			r = (r_refdef.fog_end - (-5000)) / ((5000) - (-5000));
			M_DrawSlider (220, offset+(OPT_FOG_END*8), r);

			M_Print (16, offset+(OPT_FOG_RED*8),			"        Fog red amount");
			r = (r_refdef.fog_red) / 100;
			M_DrawSlider (220, offset+(OPT_FOG_RED*8), r);

			M_Print (16, offset+(OPT_FOG_GREEN*8),			"      Fog green amount");
			r = (r_refdef.fog_green) / 100;
			M_DrawSlider (220, offset+(OPT_FOG_GREEN*8), r);

			M_Print (16, offset+(OPT_FOG_BLUE*8),			"       Fog blue amount");
			r = (r_refdef.fog_blue) / 100;
			M_DrawSlider (220, offset+(OPT_FOG_BLUE*8), r);

			M_Print (16, offset+(OPT_SCRSIZE*8),			"           Screen size");
			r = (scr_viewsize.value - 30) / (130 - 30);
			M_DrawSlider (220, offset+(OPT_SCRSIZE*8), r);

			M_Print (16, offset+(OPT_WATERTRANS*8),			"     Water tranparency");
			M_DrawSlider (220, offset+(OPT_WATERTRANS*8), r_wateralpha.value);

			M_Print (16, offset+(OPT_MIPMAP_BIAS*8),		"         MipMap Amount");
			r = (r_mipmaps_bias.value + 10) / 10;
			M_DrawSlider (220, offset+(OPT_MIPMAP_BIAS*8), r);
		
			M_Print (16, offset+(OPT_TEX_SCALEDOWN*8),		"    Texture Scale Down");
			M_DrawCheckbox (220, offset+(OPT_TEX_SCALEDOWN*8), r_tex_scale_down.value);

			M_Print (16, offset+(OPT_SMOOTH_ANIMS*8),		"Smooth Model Animation");
			M_DrawCheckbox (220, offset+(OPT_SMOOTH_ANIMS*8), r_i_model_animation.value);

			M_Print (16, offset+(OPT_SMOOTH_MOVEMENT*8),	" Smooth Model Movement");
			M_DrawCheckbox (220, offset+(OPT_SMOOTH_MOVEMENT*8), r_i_model_transform.value);

            break;

	}
	
	M_PrintWhite (16, offset+(OPT_SUBMENU*8),	 "        Select Submenu");
    switch (m_submenu)
        {
        case 0:
            M_PrintWhite (220, offset+(OPT_SUBMENU*8), "Control Options");         
            break;
        case 1:
            M_PrintWhite (220, offset+(OPT_SUBMENU*8), "Video Options");
            break;
        case 2:
            M_PrintWhite (220, offset+(OPT_SUBMENU*8), "Audio Options");
            break;
        case 3:
            M_PrintWhite (220, offset+(OPT_SUBMENU*8), "Misc. Options");
            break;
        default:
            break;
        }
}


void M_Options_Key (int k)
{
	switch (k)
	{
	case K_ESCAPE:
		M_Menu_Main_f ();
		break;

	case K_ENTER:
		m_entersound = true;
		switch (options_cursor)
		{
		case OPT_CUSTOMIZE:
			M_Menu_Keys_f ();
			break;
		case OPT_CONSOLE:
			m_state = m_none;
			Con_ToggleConsole_f ();
			break;
		case OPT_DEFAULTS1:
			Cbuf_AddText ("exec default.cfg\n");
//			Cbuf_AddText ("-klook\n");
			break;
/*
		case OPT_DEFAULTS2:
			Cbuf_AddText ("exec default2.cfg\n");
			Cbuf_AddText ("+klook\n");
			break;
		case OPT_DEFAULTS3:
			Cbuf_AddText ("exec default3.cfg\n");
			Cbuf_AddText ("+klook\n");
			break;
*/
		default:
			M_AdjustSliders (1);
			break;
		}
		return;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		options_cursor--;

		if (options_cursor == OPT_GAP)
		    options_cursor = options_cursor -1;
		    
		if (options_cursor == OPT_GAP1)
		    options_cursor = options_cursor -1;
		    
        if (m_submenu == 0)
        {
		    if (options_cursor == OPT_GAP_0)
		        options_cursor = options_cursor -1;
		    if (options_cursor == OPT_GAP_0_1)
		        options_cursor = options_cursor -1;
		    if (options_cursor == OPT_GAP_0_2)
		        options_cursor = options_cursor -1;
        }

        if (m_submenu == 1)
        {
	        if (options_cursor == OPT_GAP_1)
			    options_cursor = options_cursor -1;
	        if (options_cursor == OPT_GAP_1_1)
		    	options_cursor = options_cursor -1;
        }

        if (m_submenu == 2)
        {
	        if (options_cursor == OPT_GAP_2)
			    options_cursor = options_cursor -1;
	        if (options_cursor == OPT_GAP_2_1)
			    options_cursor = options_cursor -1;
        }

        if (m_submenu == 3)
        {
	        if (options_cursor == OPT_GAP_3)
			    options_cursor = options_cursor -1;
	        if (options_cursor == OPT_GAP_3_1)
			    options_cursor = options_cursor -1;
        }

		if (options_cursor < 0) {
			if (m_submenu == 0)
			    options_cursor = OPTIONS_ITEMS_0-1;
	        if (m_submenu == 1)
			    options_cursor = OPTIONS_ITEMS_1-1;
	        if (m_submenu == 2)
			    options_cursor = OPTIONS_ITEMS_2-1;
	        if (m_submenu == 3)
			    options_cursor = OPTIONS_ITEMS_3-1;
		}
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		options_cursor++;

		if (options_cursor == OPT_GAP)
		    options_cursor = options_cursor +1;

		if (options_cursor == OPT_GAP1)
		    options_cursor = options_cursor +1;

        if (m_submenu == 0)
        {
            if (options_cursor >= OPTIONS_ITEMS_0)
			    options_cursor = 0;
		    if (options_cursor == OPT_GAP_0)
		        options_cursor = options_cursor +1;
		    if (options_cursor == OPT_GAP_0_1)
		        options_cursor = options_cursor +1;
		    if (options_cursor == OPT_GAP_0_2)
		        options_cursor = options_cursor +1;
        }
        if (m_submenu == 1)
        {
            if (options_cursor >= OPTIONS_ITEMS_1)
			    options_cursor = 0;
	        if (options_cursor == OPT_GAP_1)
			    options_cursor = options_cursor +1;
	        if (options_cursor == OPT_GAP_1_1)
			    options_cursor = options_cursor +1;
        }
        if (m_submenu == 2)
        {
            if (options_cursor >= OPTIONS_ITEMS_2)
			    options_cursor = 0;
	        if (options_cursor == OPT_GAP_2)
			    options_cursor = options_cursor +1;
	        if (options_cursor == OPT_GAP_2_1)
			    options_cursor = options_cursor +1;
		}
        if (m_submenu == 3)
        {
            if (options_cursor >= OPTIONS_ITEMS_3)
			    options_cursor = 0;
	        if (options_cursor == OPT_GAP_3)
			    options_cursor = options_cursor +1;
	        if (options_cursor == OPT_GAP_3_1)
			    options_cursor = options_cursor +1;
		}
		break;

	case K_LEFTARROW:
		M_AdjustSliders (-1);
		break;

	case K_RIGHTARROW:
		M_AdjustSliders (1);
		break;
		
	case K_AUX1:
		m_submenu--;
		options_cursor = OPT_SUBMENU;
		if (m_submenu < 0)
		    m_submenu = 3;
		break;

	case K_AUX2:
		m_submenu++;
		options_cursor = OPT_SUBMENU;
		if (m_submenu > 3)
		    m_submenu = 0;
		break;
	}
}

//=============================================================================
/* KEYS MENU */

char *bindnames[][2] =
{
{"+attack", 		"Attack"},
{"impulse 10", 		"Next Weapon"},
{"impulse 12", 		"Previous Weapon"},
{"+jump", 			"Jump / Swim Up"},
{"+forward", 		"Move Forward"},
{"+back", 			"Move Backwards"},
{"+moveleft", 		"Move Left"},
{"+moveright", 		"Move Right"},
{"+left", 			"Turn Left"},
{"+right", 			"Turn Right"},
{"+lookup", 		"Look up"},
{"+lookdown", 		"Look down"},
{"+moveup",			"Swim Up"},
{"+movedown",		"Swim Down"},
{"+speed", 			"Run"},
{"+strafe", 		"Sidestep"},
{"centerview", 		"Center view"},
#ifdef PSP
{"+mlook", 			"Analog Nub Look"},
{"+klook", 			"D-Pad Look"},
#else
{"+mlook", 			"Mouse Look"},
{"+klook", 			"Keyboard Look"},
#endif
{"+showscores",		"Show Scores"},
{"screenshot",		"Screenshot"},
{"toggleconsole",	"Toggle Console"},
};

char *kbindnames[][2] =
{
{"+attack", 		"Attack"},
{"impulse 10", 		"Next Weapon"},
{"impulse 12", 		"Previous Weapon"},
{"impulse 13", 		"Reload / Secondary"},
{"impulse 14", 		"Zoom"},
{"+jump", 			"Jump / Swim Up"},
{"+forward", 		"Move Forward"},
{"+back", 			"Move Backwards"},
{"+moveleft", 		"Move Left"},
{"+moveright", 		"Move Right"},
{"+left", 			"Turn Left"},
{"+right", 			"Turn Right"},
{"+lookup", 		"Look up"},
{"+lookdown", 		"Look down"},
{"+moveup",			"Swim Up"},
{"+movedown",		"Swim Down"},
{"+speed", 			"Run"},
{"+strafe", 		"Sidestep"},
{"centerview", 		"Center view"},
#ifdef PSP
{"+mlook", 			"Analog Nub Look"},
{"+klook", 			"D-Pad Look"},
#else
{"+mlook", 			"Mouse Look"},
{"+klook", 			"Keyboard Look"},
#endif
{"+showscores",		"Show Scores"},
{"screenshot",		"Screenshot"},
{"toggleconsole",	"Toggle Console"},
};

#define	NUMCOMMANDS	(sizeof(bindnames)/sizeof(bindnames[0]))
#define	KNUMCOMMANDS (sizeof(kbindnames)/sizeof(kbindnames[0]))

int		keys_cursor;
int		bind_grab;

void M_Menu_Keys_f (void)
{
	key_dest = key_menu;
	m_state = m_keys;
	m_entersound = true;
}


void M_FindKeysForCommand (char *command, int *twokeys)
{
	int		count;
	int		j;
	int		l;
	char	*b;

	twokeys[0] = twokeys[1] = -1;
	l = strlen(command);
	count = 0;

	for (j=0 ; j<256 ; j++)
	{
		b = keybindings[j];
		if (!b)
			continue;
		if (!strncmp (b, command, l) )
		{
			twokeys[count] = j;
			count++;
			if (count == 2)
				break;
		}
	}
}

void M_UnbindCommand (char *command)
{
	int		j;
	int		l;
	char	*b;

	l = strlen(command);

	for (j=0 ; j<256 ; j++)
	{
		b = keybindings[j];
		if (!b)
			continue;
		if (!strncmp (b, command, l) )
			Key_SetBinding (j, "");
	}
}


void M_Keys_Draw (void)
{
	int		i, l;
	int		keys[2];
	char	*name;
	int		x, y;
	qpic_t	*p, *b;

	p = Draw_CachePic ("gfx/ttl_cstm.lmp");
	M_DrawPic ( (320-p->width)/2, 4, p);

	b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );

#ifdef PSP
	if (bind_grab)
		M_Print (12, 32, "Press a button for this action");
	else
		M_Print (18, 32, "Press CROSS to change, TRIANGLE to clear");
#else
	if (bind_grab)
		M_Print (12, 32, "Press a key or button for this action");
	else
		M_Print (18, 32, "Enter to change, backspace to clear");
#endif

// search for known bindings
    if (kurok)
    {
	    if (bind_grab)
		    M_DrawCharacter (170, 48 + keys_cursor*8, '?');
	    else
		    M_DrawCharacter (170, 48 + keys_cursor*8, 12+((int)(realtime*30)&1));
              
	for (i=0 ; i<KNUMCOMMANDS ; i++)
	{
		y = 48 + 8*i;

		M_Print (16, y, kbindnames[i][1]);

		l = strlen (kbindnames[i][0]);

		M_FindKeysForCommand (kbindnames[i][0], keys);

		if (keys[0] == -1)
		{
			M_Print (180, y, "---");
		}
		else
		{
			name = Key_KeynumToString (keys[0]);
			M_Print (180, y, name);
			x = strlen(name) * 8;
			if (keys[1] != -1)
			{
				M_Print (180 + x + 8, y, "or");
				M_Print (180 + x + 32, y, Key_KeynumToString (keys[1]));
			}
		}
	}
    }
    else
    {
	    if (bind_grab)
            M_DrawCharacter (130, 48 + keys_cursor*8, '=');
        else
            M_DrawCharacter (130, 48 + keys_cursor*8, 12+((int)(realtime*4)&1));

	for (i=0 ; i<NUMCOMMANDS ; i++)
	{
		y = 48 + 8*i;

		M_Print (16, y, bindnames[i][1]);

		l = strlen (bindnames[i][0]);

		M_FindKeysForCommand (bindnames[i][0], keys);

		if (keys[0] == -1)
		{
			M_Print (140, y, "???");
		}
		else
		{
			name = Key_KeynumToString (keys[0]);
			M_Print (140, y, name);
			x = strlen(name) * 8;
			if (keys[1] != -1)
			{
				M_Print (140 + x + 8, y, "or");
				M_Print (140 + x + 32, y, Key_KeynumToString (keys[1]));
			}
		}
	}
    }
}


void M_Keys_Key (int k)
{
	char	cmd[80];
	int		keys[2];

	if (bind_grab)
	{	// defining a key
		S_LocalSound ("misc/menu1.wav");
		if (k == K_ESCAPE)
		{
			bind_grab = false;
		}
		else if (k != '`')
		{
            if (kurok)
			    sprintf (cmd, "bind \"%s\" \"%s\"\n", Key_KeynumToString (k), kbindnames[keys_cursor][0]);
			else
			    sprintf (cmd, "bind \"%s\" \"%s\"\n", Key_KeynumToString (k), bindnames[keys_cursor][0]);
			Cbuf_InsertText (cmd);
		}

		bind_grab = false;
		return;
	}

	switch (k)
	{
	case K_ESCAPE:
		M_Menu_Options_f ();
		break;

	case K_LEFTARROW:
	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		keys_cursor--;
		if (keys_cursor < 0)
		{
		    if (kurok)
			    keys_cursor = KNUMCOMMANDS-1;
			else
			    keys_cursor = NUMCOMMANDS-1;
        }
		break;

	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound ("misc/menu1.wav");
		keys_cursor++;
		if (kurok)
		{
		    if (keys_cursor >= KNUMCOMMANDS)
			    keys_cursor = 0;
        }
        else
		{
		    if (keys_cursor >= NUMCOMMANDS)
			    keys_cursor = 0;
        }
		break;

	case K_ENTER:		// go into bind mode
	    if (kurok)
		    M_FindKeysForCommand (kbindnames[keys_cursor][0], keys);
	    else
	        M_FindKeysForCommand (bindnames[keys_cursor][0], keys);
		S_LocalSound ("misc/menu2.wav");
		if (keys[1] != -1)
		{
		    if (kurok)
			    M_UnbindCommand (kbindnames[keys_cursor][0]);
		    else
		        M_UnbindCommand (bindnames[keys_cursor][0]);
        }
		bind_grab = true;
		break;

	case K_BACKSPACE:		// delete bindings
	case K_DEL:				// delete bindings
		S_LocalSound ("misc/menu2.wav");
		if (kurok)
		    M_UnbindCommand (kbindnames[keys_cursor][0]);
        else
		    M_UnbindCommand (bindnames[keys_cursor][0]);
		break;
	}
}

//=============================================================================
/* VIDEO MENU */

void M_Menu_Video_f (void)
{
	key_dest = key_menu;
	m_state = m_video;
	m_entersound = true;
}


void M_Video_Draw (void)
{
	(*vid_menudrawfn) ();
}


void M_Video_Key (int key)
{
	(*vid_menukeyfn) (key);
}

//=============================================================================
/* HELP MENU */

int		help_page;
#define	NUM_HELP_PAGES	6
#define	KNUM_HELP_PAGES	2


void M_Menu_Help_f (void)
{
	key_dest = key_menu;
	m_state = m_help;
	m_entersound = true;
	help_page = 0;
}



void M_Help_Draw (void)
{
	if (kurok)
	    M_DrawPic (vid.width - 560, 0, Draw_CachePic ( va("gfx/menu/hp/help%i.lmp", help_page)) );
    else
        M_DrawPic (0, 0, Draw_CachePic ( va("gfx/help%i.lmp", help_page)) );
}


void M_Help_Key (int key)
{
	switch (key)
	{
	case K_ESCAPE:
		M_Menu_Main_f ();
		break;

	case K_UPARROW:
	case K_RIGHTARROW:
		m_entersound = true;
		if (!kurok)
		{
			if (++help_page >= NUM_HELP_PAGES)
				help_page = 0;
		}
		else
		{
			if (++help_page >= KNUM_HELP_PAGES)
				help_page = 0;
		}
		break;

	case K_DOWNARROW:
	case K_LEFTARROW:
		m_entersound = true;
		if (!kurok)
		{
			if (--help_page < 0)
				help_page = NUM_HELP_PAGES-1;
		}
		else
		{
			if (--help_page < 0)
				help_page = KNUM_HELP_PAGES-1;
		}
		break;
	}

}

//=============================================================================
/* QUIT MENU */

int		msgNumber;
int		m_quit_prevstate;
qboolean	wasInMenus;

#ifndef	WIN32
char *quitMessage [] = 
{
/* .........1.........2.... */
  "  Are you gonna quit    ",
  "  this game just like   ",
  "   everything else?     ",
  "                        ",
 
  " Milord, methinks that  ",
  "   thou art a lowly     ",
  " quitter. Is this true? ",
  "                        ",

  " Do I need to bust your ",
  "  face open for trying  ",
  "        to quit?        ",
  "                        ",

  " Man, I oughta smack you",
  "   for trying to quit!  ",
  "     Press Y to get     ",
  "      smacked out.      ",
 
  " Press Y to quit like a ",
  "   big loser in life.   ",
  "  Press N to stay proud ",
  "    and successful!     ",
 
  "   If you press Y to    ",
  "  quit, I will summon   ",
  "  Satan all over your   ",
  "      hard drive!       ",
 
  "  Um, Asmodeus dislikes ",
  " his children trying to ",
  " quit. Press Y to return",
  "   to your Tinkertoys.  ",
 
  "  If you quit now, I'll ",
  "  throw a blanket-party ",
  "   for you next time!   ",
  "                        "
};
#endif

void M_Menu_Quit_f (void)
{
	if (m_state == m_quit)
		return;
	wasInMenus = (key_dest == key_menu);
	key_dest = key_menu;
	m_quit_prevstate = m_state;
	m_state = m_quit;
	m_entersound = true;
	msgNumber = rand()&7;
}


void M_Quit_Key (int key)
{
	switch (key)
	{
	case K_ESCAPE:
	case 'n':
	case 'N':
		if (wasInMenus)
		{
			m_state = m_quit_prevstate;
			m_entersound = true;
		}
		else
		{
			key_dest = key_game;
			m_state = m_none;
		}
		break;

	case 'Y':
	case 'y':
#ifdef PSP
	case K_ENTER:
#endif
		key_dest = key_console;
		Host_Quit_f ();
		break;

	default:
		break;
	}

}


void M_Quit_Draw (void)
{
	if (wasInMenus)
	{
		m_state = m_quit_prevstate;
		m_recursiveDraw = true;
		M_Draw ();
		m_state = m_quit;
	}

#ifdef WIN32
	M_DrawTextBox (0, 0, 38, 23);
	M_PrintWhite (16, 12,  "  Quake version 1.09 by id Software\n\n");
	M_PrintWhite (16, 28,  "Programming        Art \n");
	M_Print (16, 36,  " John Carmack       Adrian Carmack\n");
	M_Print (16, 44,  " Michael Abrash     Kevin Cloud\n");
	M_Print (16, 52,  " John Cash          Paul Steed\n");
	M_Print (16, 60,  " Dave 'Zoid' Kirsch\n");
	M_PrintWhite (16, 68,  "Design             Biz\n");
	M_Print (16, 76,  " John Romero        Jay Wilbur\n");
	M_Print (16, 84,  " Sandy Petersen     Mike Wilson\n");
	M_Print (16, 92,  " American McGee     Donna Jackson\n");
	M_Print (16, 100,  " Tim Willits        Todd Hollenshead\n");
	M_PrintWhite (16, 108, "Support            Projects\n");
	M_Print (16, 116, " Barrett Alexander  Shawn Green\n");
	M_PrintWhite (16, 124, "Sound Effects\n");
	M_Print (16, 132, " Trent Reznor and Nine Inch Nails\n\n");
	M_PrintWhite (16, 140, "Quake is a trademark of Id Software,\n");
	M_PrintWhite (16, 148, "inc., (c)1996 Id Software, inc. All\n");
	M_PrintWhite (16, 156, "rights reserved. NIN logo is a\n");
	M_PrintWhite (16, 164, "registered trademark licensed to\n");
	M_PrintWhite (16, 172, "Nothing Interactive, Inc. All rights\n");
	M_PrintWhite (16, 180, "reserved. Press y to exit\n");
#elif defined PSP
	M_DrawTextBox (56, 76, 24, 4);
	M_Print (64, 84,	"      Really quit?      ");
	M_Print (64, 92,	"                        ");
	M_Print (64, 100,	"  Press CROSS to quit,  ");
	M_Print (64, 108,	" or CIRCLE to continue. ");
#else
	M_DrawTextBox (56, 76, 24, 4);
	M_Print (64, 84,  quitMessage[msgNumber*4+0]);
	M_Print (64, 92,  quitMessage[msgNumber*4+1]);
	M_Print (64, 100, quitMessage[msgNumber*4+2]);
	M_Print (64, 108, quitMessage[msgNumber*4+3]);
#endif
}
//=============================================================================
/* OSK IMPLEMENTATION */
#define CHAR_SIZE 8
#define MAX_Y 8
#define MAX_X 12

#define MAX_CHAR_LINE 36
#define MAX_CHAR      72

int  osk_pos_x = 0;
int  osk_pos_y = 0;
int  max_len   = 0;
int  m_old_state = 0;

char* osk_out_buff = NULL;
char  osk_buffer[128];

char *osk_text [] = 
	{ 
		" 1 2 3 4 5 6 7 8 9 0 - = ` ",
		" q w e r t y u i o p [ ]   ",
		"   a s d f g h j k l ; ' \\ ",
		"     z x c v b n m   , . / ",
		"                           ",
		" ! @ # $ % ^ & * ( ) _ + ~ ",
		" Q W E R T Y U I O P { }   ",
		"   A S D F G H J K L : \" | ",
		"     Z X C V B N M   < > ? "
	};

char *osk_help [] = 
	{ 
		"CONFIRM: ",
		" SQUARE  ",
		"CANCEL:  ",
		" CIRCLE  ",
		"DELETE:  ",
		" TRIAGLE ",
		"ADD CHAR:",
		" CROSS   ",
		""
	};

void M_Menu_OSK_f (char *input, char *output, int outlen)
{
	key_dest = key_menu;
	m_old_state = m_state;
	m_state = m_osk;
	m_entersound = false;
	max_len = outlen;
	strncpy(osk_buffer,input,max_len);
	osk_buffer[outlen] = '\0';
	osk_out_buff = output; 
}

void Con_OSK_f (char *input, char *output, int outlen)
{
	max_len = outlen;
	strncpy(osk_buffer,input,max_len);
	osk_buffer[outlen] = '\0';
	osk_out_buff = output; 
}


void M_OSK_Draw (void)
{
#ifdef PSP
	int x,y;
	int i;
	
	char *selected_line = osk_text[osk_pos_y]; 
	char selected_char[2];
	
	selected_char[0] = selected_line[1+(2*osk_pos_x)];
	selected_char[1] = '\0';
	if (selected_char[0] == ' ' || selected_char[0] == '\t') 
		selected_char[0] = 'X';
		
	y = 20;
	x = 16;

	M_DrawTextBox (10, 10, 		     26, 10);
	M_DrawTextBox (10+(26*CHAR_SIZE),    10,  10, 10);
	M_DrawTextBox (10, 10+(10*CHAR_SIZE),36,  3);
	
	for(i=0;i<=MAX_Y;i++) 
	{
		M_PrintWhite (x, y+(CHAR_SIZE*i), osk_text[i]);
		if (i % 2 == 0)
			M_Print      (x+(27*CHAR_SIZE), y+(CHAR_SIZE*i), osk_help[i]);
		else			
			M_PrintWhite (x+(27*CHAR_SIZE), y+(CHAR_SIZE*i), osk_help[i]);
	}
	
	int text_len = strlen(osk_buffer);
	if (text_len > MAX_CHAR_LINE) {
		
		char oneline[MAX_CHAR_LINE+1];
		strncpy(oneline,osk_buffer,MAX_CHAR_LINE);
		oneline[MAX_CHAR_LINE] = '\0';
		
		M_Print (x+4, y+4+(CHAR_SIZE*(MAX_Y+2)), oneline );
		
		strncpy(oneline,osk_buffer+MAX_CHAR_LINE, text_len - MAX_CHAR_LINE);
		oneline[text_len - MAX_CHAR_LINE] = '\0';
		
		M_Print (x+4, y+4+(CHAR_SIZE*(MAX_Y+3)), oneline );
		M_PrintWhite (x+4+(CHAR_SIZE*(text_len - MAX_CHAR_LINE)), y+4+(CHAR_SIZE*(MAX_Y+3)),"_");
	}
	else {
		M_Print (x+4, y+4+(CHAR_SIZE*(MAX_Y+2)), osk_buffer );
		M_PrintWhite (x+4+(CHAR_SIZE*(text_len)), y+4+(CHAR_SIZE*(MAX_Y+2)),"_");
	}
	M_Print      (x+((((osk_pos_x)*2)+1)*CHAR_SIZE), y+(osk_pos_y*CHAR_SIZE), selected_char);

#endif
}

void M_OSK_Key (int key)
{
#ifdef PSP
	switch (key)
	{
	case K_RIGHTARROW:
		osk_pos_x++;
		if (osk_pos_x > MAX_X)
			osk_pos_x = MAX_X;
		break;
	case K_LEFTARROW:
		osk_pos_x--;
		if (osk_pos_x < 0)
			osk_pos_x = 0;
		break;
	case K_DOWNARROW:
		osk_pos_y++;
		if (osk_pos_y > MAX_Y)
			osk_pos_y = MAX_Y;
		break;
	case K_UPARROW:
		osk_pos_y--;
		if (osk_pos_y < 0)
			osk_pos_y = 0;
		break;
	case K_ENTER: 
		if (max_len > strlen(osk_buffer)) {
			char *selected_line = osk_text[osk_pos_y]; 
			char selected_char[2];
			
			selected_char[0] = selected_line[1+(2*osk_pos_x)];
			
			if (selected_char[0] == '\t')
				selected_char[0] = ' ';
			
			selected_char[1] = '\0';
			strcat(osk_buffer,selected_char);		
		}
		break;
	case K_DEL:
		if (strlen(osk_buffer) > 0) {
			osk_buffer[strlen(osk_buffer)-1] = '\0';	
		}
		break;
	case K_INS:
		strncpy(osk_out_buff,osk_buffer,max_len);
		
		m_state = m_old_state;
		break;
	case K_ESCAPE:
		m_state = m_old_state;
		break;
	default:
		break;
	}
#endif		
}

void Con_OSK_Key (int key)
{
#ifdef PSP
	switch (key)
	{
	case K_RIGHTARROW:
		osk_pos_x++;
		if (osk_pos_x > MAX_X)
			osk_pos_x = MAX_X;
		break;
	case K_LEFTARROW:
		osk_pos_x--;
		if (osk_pos_x < 0)
			osk_pos_x = 0;
		break;
	case K_DOWNARROW:
		osk_pos_y++;
		if (osk_pos_y > MAX_Y)
			osk_pos_y = MAX_Y;
		break;
	case K_UPARROW:
		osk_pos_y--;
		if (osk_pos_y < 0)
			osk_pos_y = 0;
		break;
	case K_ENTER: 
		if (max_len > strlen(osk_buffer)) {
			char *selected_line = osk_text[osk_pos_y]; 
			char selected_char[2];
			
			selected_char[0] = selected_line[1+(2*osk_pos_x)];
			
			if (selected_char[0] == '\t')
				selected_char[0] = ' ';
			
			selected_char[1] = '\0';
			strcat(osk_buffer,selected_char);		
		}
		break;
	case K_DEL:
		if (strlen(osk_buffer) > 0) {
			osk_buffer[strlen(osk_buffer)-1] = '\0';	
		}
		break;
	case K_INS:
		strncpy(osk_out_buff,osk_buffer,max_len);
		Con_SetOSKActive(false);
		break;
	case K_ESCAPE:
		Con_SetOSKActive(false);
		break;
	default:
		break;
	}
#endif		
}	

//=============================================================================

/* SERIAL CONFIG MENU */

int		serialConfig_cursor;
int		serialConfig_cursor_table[] = {48, 64, 80, 96, 112, 132};
#define	NUM_SERIALCONFIG_CMDS	6

static int ISA_uarts[]	= {0x3f8,0x2f8,0x3e8,0x2e8};
static int ISA_IRQs[]	= {4,3,4,3};
int serialConfig_baudrate[] = {9600,14400,19200,28800,38400,57600};

int		serialConfig_comport;
int		serialConfig_irq ;
int		serialConfig_baud;
char	serialConfig_phone[16];

void M_Menu_SerialConfig_f (void)
{
	int		n;
	int		port;
	int		baudrate;
	qboolean	useModem;

	key_dest = key_menu;
	m_state = m_serialconfig;
	m_entersound = true;
	if (JoiningGame && SerialConfig)
		serialConfig_cursor = 4;
	else
		serialConfig_cursor = 5;

	(*GetComPortConfig) (0, &port, &serialConfig_irq, &baudrate, &useModem);

	// map uart's port to COMx
	for (n = 0; n < 4; n++)
		if (ISA_uarts[n] == port)
			break;
	if (n == 4)
	{
		n = 0;
		serialConfig_irq = 4;
	}
	serialConfig_comport = n + 1;

	// map baudrate to index
	for (n = 0; n < 6; n++)
		if (serialConfig_baudrate[n] == baudrate)
			break;
	if (n == 6)
		n = 5;
	serialConfig_baud = n;

	m_return_onerror = false;
	m_return_reason[0] = 0;
}


void M_SerialConfig_Draw (void)
{
	qpic_t	*p,*b;
	int		basex;
	char	*startJoin;
	char	*directModem;
	
	if (!kurok)
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );

    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
	p = Draw_CachePic ("gfx/p_multi.lmp");
	basex = (320-p->width)/2;
	M_DrawPic (basex, 4, p);

	if (StartingGame)
		startJoin = "New Game";
	else
		startJoin = "Join Game";
	if (SerialConfig)
		directModem = "Modem";
	else
		directModem = "Direct Connect";
	M_Print (basex, 32, va ("%s - %s", startJoin, directModem));
	basex += 8;

	M_Print (basex, serialConfig_cursor_table[0], "Port");
	M_DrawTextBox (160, 40, 4, 1);
	M_Print (168, serialConfig_cursor_table[0], va("COM%u", serialConfig_comport));

	M_Print (basex, serialConfig_cursor_table[1], "IRQ");
	M_DrawTextBox (160, serialConfig_cursor_table[1]-8, 1, 1);
	M_Print (168, serialConfig_cursor_table[1], va("%u", serialConfig_irq));

	M_Print (basex, serialConfig_cursor_table[2], "Baud");
	M_DrawTextBox (160, serialConfig_cursor_table[2]-8, 5, 1);
	M_Print (168, serialConfig_cursor_table[2], va("%u", serialConfig_baudrate[serialConfig_baud]));

	if (SerialConfig)
	{
		M_Print (basex, serialConfig_cursor_table[3], "Modem Setup...");
		if (JoiningGame)
		{
			M_Print (basex, serialConfig_cursor_table[4], "Phone number");
			M_DrawTextBox (160, serialConfig_cursor_table[4]-8, 16, 1);
			M_Print (168, serialConfig_cursor_table[4], serialConfig_phone);
		}
	}

	if (JoiningGame)
	{
		M_DrawTextBox (basex, serialConfig_cursor_table[5]-8, 7, 1);
		M_Print (basex+8, serialConfig_cursor_table[5], "Connect");
	}
	else
	{
		M_DrawTextBox (basex, serialConfig_cursor_table[5]-8, 2, 1);
		M_Print (basex+8, serialConfig_cursor_table[5], "OK");
	}

	M_DrawCharacter (basex-8, serialConfig_cursor_table [serialConfig_cursor], 12+((int)(realtime*4)&1));

	if (serialConfig_cursor == 4)
		M_DrawCharacter (168 + 8*strlen(serialConfig_phone), serialConfig_cursor_table [serialConfig_cursor], 10+((int)(realtime*4)&1));

	if (*m_return_reason)
		M_PrintWhite (basex, 148, m_return_reason);
}


void M_SerialConfig_Key (int key)
{
	int		l;

	switch (key)
	{
	case K_ESCAPE:
		M_Menu_Net_f ();
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		serialConfig_cursor--;
		if (serialConfig_cursor < 0)
			serialConfig_cursor = NUM_SERIALCONFIG_CMDS-1;
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		serialConfig_cursor++;
		if (serialConfig_cursor >= NUM_SERIALCONFIG_CMDS)
			serialConfig_cursor = 0;
		break;

	case K_LEFTARROW:
		if (serialConfig_cursor > 2)
			break;
		S_LocalSound ("misc/menu3.wav");

		if (serialConfig_cursor == 0)
		{
			serialConfig_comport--;
			if (serialConfig_comport == 0)
				serialConfig_comport = 4;
			serialConfig_irq = ISA_IRQs[serialConfig_comport-1];
		}

		if (serialConfig_cursor == 1)
		{
			serialConfig_irq--;
			if (serialConfig_irq == 6)
				serialConfig_irq = 5;
			if (serialConfig_irq == 1)
				serialConfig_irq = 7;
		}

		if (serialConfig_cursor == 2)
		{
			serialConfig_baud--;
			if (serialConfig_baud < 0)
				serialConfig_baud = 5;
		}

		break;

	case K_RIGHTARROW:
		if (serialConfig_cursor > 2)
			break;
forward:
		S_LocalSound ("misc/menu3.wav");

		if (serialConfig_cursor == 0)
		{
			serialConfig_comport++;
			if (serialConfig_comport > 4)
				serialConfig_comport = 1;
			serialConfig_irq = ISA_IRQs[serialConfig_comport-1];
		}

		if (serialConfig_cursor == 1)
		{
			serialConfig_irq++;
			if (serialConfig_irq == 6)
				serialConfig_irq = 7;
			if (serialConfig_irq == 8)
				serialConfig_irq = 2;
		}

		if (serialConfig_cursor == 2)
		{
			serialConfig_baud++;
			if (serialConfig_baud > 5)
				serialConfig_baud = 0;
		}

		break;

	case K_ENTER:
		if (serialConfig_cursor < 3)
			goto forward;

		m_entersound = true;

		if (serialConfig_cursor == 3)
		{
			(*SetComPortConfig) (0, ISA_uarts[serialConfig_comport-1], serialConfig_irq, serialConfig_baudrate[serialConfig_baud], SerialConfig);

			M_Menu_ModemConfig_f ();
			break;
		}

		if (serialConfig_cursor == 4)
		{
			serialConfig_cursor = 5;
			break;
		}

		// serialConfig_cursor == 5 (OK/CONNECT)
		(*SetComPortConfig) (0, ISA_uarts[serialConfig_comport-1], serialConfig_irq, serialConfig_baudrate[serialConfig_baud], SerialConfig);

		M_ConfigureNetSubsystem ();

		if (StartingGame)
		{
			M_Menu_GameOptions_f ();
			break;
		}

		m_return_state = m_state;
		m_return_onerror = true;
		key_dest = key_game;
		m_state = m_none;

		if (SerialConfig)
			Cbuf_AddText (va ("connect \"%s\"\n", serialConfig_phone));
		else
			Cbuf_AddText ("connect\n");
		break;

	case K_BACKSPACE:
		if (serialConfig_cursor == 4)
		{
			if (strlen(serialConfig_phone))
				serialConfig_phone[strlen(serialConfig_phone)-1] = 0;
		}
		break;

	default:
		if (key < 32 || key > 127)
			break;
		if (serialConfig_cursor == 4)
		{
			l = strlen(serialConfig_phone);
			if (l < 15)
			{
				serialConfig_phone[l+1] = 0;
				serialConfig_phone[l] = key;
			}
		}
	}

	if (DirectConfig && (serialConfig_cursor == 3 || serialConfig_cursor == 4))
	{
		if (key == K_UPARROW)
			serialConfig_cursor = 2;
		else
			serialConfig_cursor = 5;
	}
	if (SerialConfig && StartingGame && serialConfig_cursor == 4)
	{
		if (key == K_UPARROW)
			serialConfig_cursor = 3;
		else
			serialConfig_cursor = 5;
	}
}

//=============================================================================
/* MODEM CONFIG MENU */

int		modemConfig_cursor;
int		modemConfig_cursor_table [] = {40, 56, 88, 120, 156};
#define NUM_MODEMCONFIG_CMDS	5

char	modemConfig_dialing;
char	modemConfig_clear [16];
char	modemConfig_init [32];
char	modemConfig_hangup [16];

void M_Menu_ModemConfig_f (void)
{
	key_dest = key_menu;
	m_state = m_modemconfig;
	m_entersound = true;
	(*GetModemConfig) (0, &modemConfig_dialing, modemConfig_clear, modemConfig_init, modemConfig_hangup);
}


void M_ModemConfig_Draw (void)
{
	qpic_t	*p,*b;
	int		basex;
	
	if (!kurok)
		M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );

    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
	p = Draw_CachePic ("gfx/p_multi.lmp");
	basex = (320-p->width)/2;
	M_DrawPic (basex, 4, p);
	basex += 8;

	if (modemConfig_dialing == 'P')
		M_Print (basex, modemConfig_cursor_table[0], "Pulse Dialing");
	else
		M_Print (basex, modemConfig_cursor_table[0], "Touch Tone Dialing");

	M_Print (basex, modemConfig_cursor_table[1], "Clear");
	M_DrawTextBox (basex, modemConfig_cursor_table[1]+4, 16, 1);
	M_Print (basex+8, modemConfig_cursor_table[1]+12, modemConfig_clear);
	if (modemConfig_cursor == 1)
		M_DrawCharacter (basex+8 + 8*strlen(modemConfig_clear), modemConfig_cursor_table[1]+12, 10+((int)(realtime*4)&1));

	M_Print (basex, modemConfig_cursor_table[2], "Init");
	M_DrawTextBox (basex, modemConfig_cursor_table[2]+4, 30, 1);
	M_Print (basex+8, modemConfig_cursor_table[2]+12, modemConfig_init);
	if (modemConfig_cursor == 2)
		M_DrawCharacter (basex+8 + 8*strlen(modemConfig_init), modemConfig_cursor_table[2]+12, 10+((int)(realtime*4)&1));

	M_Print (basex, modemConfig_cursor_table[3], "Hangup");
	M_DrawTextBox (basex, modemConfig_cursor_table[3]+4, 16, 1);
	M_Print (basex+8, modemConfig_cursor_table[3]+12, modemConfig_hangup);
	if (modemConfig_cursor == 3)
		M_DrawCharacter (basex+8 + 8*strlen(modemConfig_hangup), modemConfig_cursor_table[3]+12, 10+((int)(realtime*4)&1));

	M_DrawTextBox (basex, modemConfig_cursor_table[4]-8, 2, 1);
	M_Print (basex+8, modemConfig_cursor_table[4], "OK");

	M_DrawCharacter (basex-8, modemConfig_cursor_table [modemConfig_cursor], 12+((int)(realtime*4)&1));
}


void M_ModemConfig_Key (int key)
{
	int		l;

	switch (key)
	{
	case K_ESCAPE:
		M_Menu_SerialConfig_f ();
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		modemConfig_cursor--;
		if (modemConfig_cursor < 0)
			modemConfig_cursor = NUM_MODEMCONFIG_CMDS-1;
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		modemConfig_cursor++;
		if (modemConfig_cursor >= NUM_MODEMCONFIG_CMDS)
			modemConfig_cursor = 0;
		break;

	case K_LEFTARROW:
	case K_RIGHTARROW:
		if (modemConfig_cursor == 0)
		{
			if (modemConfig_dialing == 'P')
				modemConfig_dialing = 'T';
			else
				modemConfig_dialing = 'P';
			S_LocalSound ("misc/menu1.wav");
		}
		break;

	case K_ENTER:
		if (modemConfig_cursor == 0)
		{
			if (modemConfig_dialing == 'P')
				modemConfig_dialing = 'T';
			else
				modemConfig_dialing = 'P';
			m_entersound = true;
		}

		if (modemConfig_cursor == 4)
		{
			(*SetModemConfig) (0, va ("%c", modemConfig_dialing), modemConfig_clear, modemConfig_init, modemConfig_hangup);
			m_entersound = true;
			M_Menu_SerialConfig_f ();
		}
		break;

	case K_BACKSPACE:
		if (modemConfig_cursor == 1)
		{
			if (strlen(modemConfig_clear))
				modemConfig_clear[strlen(modemConfig_clear)-1] = 0;
		}

		if (modemConfig_cursor == 2)
		{
			if (strlen(modemConfig_init))
				modemConfig_init[strlen(modemConfig_init)-1] = 0;
		}

		if (modemConfig_cursor == 3)
		{
			if (strlen(modemConfig_hangup))
				modemConfig_hangup[strlen(modemConfig_hangup)-1] = 0;
		}
		break;

	default:
		if (key < 32 || key > 127)
			break;

		if (modemConfig_cursor == 1)
		{
			l = strlen(modemConfig_clear);
			if (l < 15)
			{
				modemConfig_clear[l+1] = 0;
				modemConfig_clear[l] = key;
			}
		}

		if (modemConfig_cursor == 2)
		{
			l = strlen(modemConfig_init);
			if (l < 29)
			{
				modemConfig_init[l+1] = 0;
				modemConfig_init[l] = key;
			}
		}

		if (modemConfig_cursor == 3)
		{
			l = strlen(modemConfig_hangup);
			if (l < 15)
			{
				modemConfig_hangup[l+1] = 0;
				modemConfig_hangup[l] = key;
			}
		}
	}
}

//=============================================================================
/* LAN CONFIG MENU */

int		lanConfig_cursor = -1;
int		lanConfig_cursor_table [] = {72, 92, 124};
#define NUM_LANCONFIG_CMDS	3

int 	lanConfig_port;
char	lanConfig_portname[6];
char	lanConfig_joinname[22];

void M_Menu_LanConfig_f (void)
{
	key_dest = key_menu;
	m_state = m_lanconfig;
	m_entersound = true;
	if (lanConfig_cursor == -1)
	{
		if (JoiningGame && TCPIPConfig)
			lanConfig_cursor = 2;
		else
			lanConfig_cursor = 1;
	}
	if (StartingGame && lanConfig_cursor == 2)
		lanConfig_cursor = 1;
	lanConfig_port = DEFAULTnet_hostport;
	sprintf(lanConfig_portname, "%u", lanConfig_port);

	m_return_onerror = false;
	m_return_reason[0] = 0;
}


void M_LanConfig_Draw (void)
{
	qpic_t	*p,*b;
	int		basex;
	char	*startJoin;
	char	*protocol;
	
	if (!kurok)
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );

    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
	p = Draw_CachePic ("gfx/p_multi.lmp");
	basex = (320-p->width)/2;
	M_DrawPic (basex, 4, p);

	if (StartingGame)
		startJoin = "New Game";
	else
		startJoin = "Join Game";
	if (IPXConfig)
		protocol = "IPX";
	else
		protocol = "TCP/IP";
	M_PrintWhite (basex, 32, va ("%s - %s", startJoin, protocol));
	basex += 8;

	M_PrintWhite (basex, 52, "Address:");
	if (IPXConfig)
		M_Print (basex+9*8, 52, my_ipx_address);
	else
		M_Print (basex+9*8, 52, my_tcpip_address);

	M_PrintWhite (basex, lanConfig_cursor_table[0], "Port");
	M_DrawTextBox (basex+8*8, lanConfig_cursor_table[0]-8, 6, 1);
	M_PrintWhite (basex+9*8, lanConfig_cursor_table[0], lanConfig_portname);

	if (JoiningGame)
	{
		M_PrintWhite (basex, lanConfig_cursor_table[1], "Search for local games...");
		M_PrintWhite (basex, 108, "Join game at:");
		M_DrawTextBox (basex+8, lanConfig_cursor_table[2]-8, 22, 1);
		M_PrintWhite (basex+16, lanConfig_cursor_table[2], lanConfig_joinname);
	}
	else
	{
		M_DrawTextBox (basex, lanConfig_cursor_table[1]-8, 2, 1);
		M_PrintWhite (basex+8, lanConfig_cursor_table[1], "OK");
	}

	if (*m_return_reason)
		M_PrintWhite (basex, 148, m_return_reason);

	if (kurok)
	{
	    M_DrawCharacter (basex-8, lanConfig_cursor_table [lanConfig_cursor], 12+((int)(realtime*30)&1));

	    if (lanConfig_cursor == 0)
		    M_DrawCharacter (basex+9*8 + 8*strlen(lanConfig_portname), lanConfig_cursor_table [0], 10+((int)(realtime*30)&1));

	    if (lanConfig_cursor == 2)
		    M_DrawCharacter (basex+16 + 8*strlen(lanConfig_joinname), lanConfig_cursor_table [2], 10+((int)(realtime*30)&1));
    }
	else
	{
	    M_DrawCharacter (basex-8, lanConfig_cursor_table [lanConfig_cursor], 12+((int)(realtime*4)&1));

	    if (lanConfig_cursor == 0)
		    M_DrawCharacter (basex+9*8 + 8*strlen(lanConfig_portname), lanConfig_cursor_table [0], 10+((int)(realtime*4)&1));

	    if (lanConfig_cursor == 2)
		    M_DrawCharacter (basex+16 + 8*strlen(lanConfig_joinname), lanConfig_cursor_table [2], 10+((int)(realtime*4)&1));
    }

}


void M_LanConfig_Key (int key)
{
	int		l;

	switch (key)
	{
	case K_ESCAPE:
//		M_Menu_Net_f ();
		M_Menu_MultiPlayer_f ();
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		lanConfig_cursor--;
		if (lanConfig_cursor < 0)
			lanConfig_cursor = NUM_LANCONFIG_CMDS-1;
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		lanConfig_cursor++;
		if (lanConfig_cursor >= NUM_LANCONFIG_CMDS)
			lanConfig_cursor = 0;
		break;

	case K_INS:
		if (lanConfig_cursor == 0)
		{
			M_Menu_OSK_f(lanConfig_portname, lanConfig_portname, 6);
			break;
		}

		if (lanConfig_cursor == 2)
		{
			M_Menu_OSK_f(lanConfig_joinname, lanConfig_joinname, 22);
			break;
		}
		break;

	case K_ENTER:
		if (lanConfig_cursor == 0)
			break;

		m_entersound = true;

		M_ConfigureNetSubsystem ();

		if (lanConfig_cursor == 1)
		{
			if (StartingGame)
			{
				M_Menu_GameOptions_f ();
				break;
			}
			M_Menu_Search_f();
			break;
		}

		if (lanConfig_cursor == 2)
		{
			m_return_state = m_state;
			m_return_onerror = true;
			key_dest = key_game;
			m_state = m_none;
			Cbuf_AddText ( va ("connect \"%s\"\n", lanConfig_joinname) );
			break;
		}

		break;

	case K_BACKSPACE:
		if (lanConfig_cursor == 0)
		{
			if (strlen(lanConfig_portname))
				lanConfig_portname[strlen(lanConfig_portname)-1] = 0;
		}

		if (lanConfig_cursor == 2)
		{
			if (strlen(lanConfig_joinname))
				lanConfig_joinname[strlen(lanConfig_joinname)-1] = 0;
		}
		break;

	default:
		if (key < 32 || key > 127)
			break;

		if (lanConfig_cursor == 2)
		{
			l = strlen(lanConfig_joinname);
			if (l < 21)
			{
				lanConfig_joinname[l+1] = 0;
				lanConfig_joinname[l] = key;
			}
		}

		if (key < '0' || key > '9')
			break;
		if (lanConfig_cursor == 0)
		{
			l = strlen(lanConfig_portname);
			if (l < 5)
			{
				lanConfig_portname[l+1] = 0;
				lanConfig_portname[l] = key;
			}
		}
	}

	if (StartingGame && lanConfig_cursor == 2)
	{
		if (key == K_UPARROW)
			lanConfig_cursor = 1;
		else
			lanConfig_cursor = 0;
	}

	l =  Q_atoi(lanConfig_portname);
	if (l > 65535)
		l = lanConfig_port;
	else
		lanConfig_port = l;
	sprintf(lanConfig_portname, "%u", lanConfig_port);
}

//=============================================================================
/* GAME OPTIONS MENU */

typedef struct
{
	char	*name;
	char	*description;
} level_t;

level_t		levels[] =
{
	{"start", "Entrance"},	// 0

	{"e1m1", "Slipgate Complex"},				// 1
	{"e1m2", "Castle of the Damned"},
	{"e1m3", "The Necropolis"},
	{"e1m4", "The Grisly Grotto"},
	{"e1m5", "Gloom Keep"},
	{"e1m6", "The Door To Chthon"},
	{"e1m7", "The House of Chthon"},
	{"e1m8", "Ziggurat Vertigo"},

	{"e2m1", "The Installation"},				// 9
	{"e2m2", "Ogre Citadel"},
	{"e2m3", "Crypt of Decay"},
	{"e2m4", "The Ebon Fortress"},
	{"e2m5", "The Wizard's Manse"},
	{"e2m6", "The Dismal Oubliette"},
	{"e2m7", "Underearth"},

	{"e3m1", "Termination Central"},			// 16
	{"e3m2", "The Vaults of Zin"},
	{"e3m3", "The Tomb of Terror"},
	{"e3m4", "Satan's Dark Delight"},
	{"e3m5", "Wind Tunnels"},
	{"e3m6", "Chambers of Torment"},
	{"e3m7", "The Haunted Halls"},

	{"e4m1", "The Sewage System"},				// 23
	{"e4m2", "The Tower of Despair"},
	{"e4m3", "The Elder God Shrine"},
	{"e4m4", "The Palace of Hate"},
	{"e4m5", "Hell's Atrium"},
	{"e4m6", "The Pain Maze"},
	{"e4m7", "Azure Agony"},
	{"e4m8", "The Nameless City"},

	{"end", "Shub-Niggurath's Pit"},			// 31

	{"dm1", "Place of Two Deaths"},				// 32
	{"dm2", "Claustrophobopolis"},
	{"dm3", "The Abandoned Base"},
	{"dm4", "The Bad Place"},
	{"dm5", "The Cistern"},
	{"dm6", "The Dark Zone"}
};

level_t		kuroklevels[] =
{
	{"start", "Entrance"},		// 0
	{"e1m1", "Base Entrance"},
	{"e1m2", "Base"},
	{"e1m3", "Canyon Testing Grounds"},
	{"e1m4", "Cavern Testing Grounds"},
	{"e1m5", "Underground Base"},
	{"e1m6", "Experiment Rex"},

	{"kdm1", "Canyon Arena"}, 	// 7
	{"kdm2", "Base Arena"},
	{"kdm3", "Ruins Arena"},
	{"kdm4", "Classic Complex"}
};

//MED 01/06/97 added hipnotic levels
level_t     hipnoticlevels[] =
{
   {"start", "Command HQ"},  // 0

   {"hip1m1", "The Pumping Station"},          // 1
   {"hip1m2", "Storage Facility"},
   {"hip1m3", "The Lost Mine"},
   {"hip1m4", "Research Facility"},
   {"hip1m5", "Military Complex"},

   {"hip2m1", "Ancient Realms"},          // 6
   {"hip2m2", "The Black Cathedral"},
   {"hip2m3", "The Catacombs"},
   {"hip2m4", "The Crypt"},
   {"hip2m5", "Mortum's Keep"},
   {"hip2m6", "The Gremlin's Domain"},

   {"hip3m1", "Tur Torment"},       // 12
   {"hip3m2", "Pandemonium"},
   {"hip3m3", "Limbo"},
   {"hip3m4", "The Gauntlet"},

   {"hipend", "Armagon's Lair"},       // 16

   {"hipdm1", "The Edge of Oblivion"}           // 17
};

//PGM 01/07/97 added rogue levels
//PGM 03/02/97 added dmatch level
level_t		roguelevels[] =
{
	{"start",	"Split Decision"},
	{"r1m1",	"Deviant's Domain"},
	{"r1m2",	"Dread Portal"},
	{"r1m3",	"Judgement Call"},
	{"r1m4",	"Cave of Death"},
	{"r1m5",	"Towers of Wrath"},
	{"r1m6",	"Temple of Pain"},
	{"r1m7",	"Tomb of the Overlord"},
	{"r2m1",	"Tempus Fugit"},
	{"r2m2",	"Elemental Fury I"},
	{"r2m3",	"Elemental Fury II"},
	{"r2m4",	"Curse of Osiris"},
	{"r2m5",	"Wizard's Keep"},
	{"r2m6",	"Blood Sacrifice"},
	{"r2m7",	"Last Bastion"},
	{"r2m8",	"Source of Evil"},
	{"ctf1",    "Division of Change"}
};

typedef struct
{
	char	*description;
	int		firstLevel;
	int		levels;
} episode_t;

episode_t	episodes[] =
{
	{"Welcome to Quake", 0, 1},
	{"Doomed Dimension", 1, 8},
	{"Realm of Black Magic", 9, 7},
	{"Netherworld", 16, 7},
	{"The Elder World", 23, 8},
	{"Final Level", 31, 1},
	{"Deathmatch Arena", 32, 6}
};

episode_t	kurokepisodes[] =
{
	{"Kurok Hub", 0, 1},
	{"Jungle Base Chapter", 1, 6},
	{"Kurok Arena", 7, 4}
};

//MED 01/06/97  added hipnotic episodes
episode_t   hipnoticepisodes[] =
{
   {"Scourge of Armagon", 0, 1},
   {"Fortress of the Dead", 1, 5},
   {"Dominion of Darkness", 6, 6},
   {"The Rift", 12, 4},
   {"Final Level", 16, 1},
   {"Deathmatch Arena", 17, 1}
};

//PGM 01/07/97 added rogue episodes
//PGM 03/02/97 added dmatch episode
episode_t	rogueepisodes[] =
{
	{"Introduction", 0, 1},
	{"Hell's Fortress", 1, 7},
	{"Corridors of Time", 8, 8},
	{"Deathmatch Arena", 16, 1}
};

int	startepisode;
int	startlevel;
int maxplayers;
qboolean m_serverInfoMessage = false;
double m_serverInfoMessageTime;

void M_Menu_GameOptions_f (void)
{
	key_dest = key_menu;
	m_state = m_gameoptions;
	m_entersound = true;
	if (maxplayers == 0)
		maxplayers = svs.maxclients;
	if (maxplayers < 2)
		maxplayers = svs.maxclientslimit;
}


int gameoptions_cursor_table[] = {40, 56, 64, 72, 80, 88, 96, 104, 112, 128, 136};
#define	NUM_GAMEOPTIONS	11
int		gameoptions_cursor;

void M_GameOptions_Draw (void)
{
	qpic_t	*p,*b;
	int		x;
//	int		r;
	
	if (kurok)
        // line cursor
	    M_DrawCharacter (144, gameoptions_cursor_table[gameoptions_cursor], 12+((int)(realtime*30)&1));
	else
	{
	    M_DrawTransPic (16, 4, Draw_CachePic ("gfx/qplaque.lmp") );
        // line cursor
	    M_DrawCharacter (144, gameoptions_cursor_table[gameoptions_cursor], 12+((int)(realtime*4)&1));
    }

    b = Draw_CachePic ("gfx/m_bttns.lmp");
	M_DrawPic ( (320-b->width)/2, 248, b );
	p = Draw_CachePic ("gfx/p_multi.lmp");
	M_DrawPic ( (320-p->width)/2, 4, p);

	M_DrawTextBox (152, 32, 10, 1);
	M_PrintWhite (160, 40, "begin game");

	M_PrintWhite (0, 56, "      Max players");
	M_Print (160, 56, va("%i", maxplayers) );

	M_PrintWhite (0, 64, "        Game Type");
	if (coop.value)
		M_Print (160, 64, "Cooperative");
	else
		M_Print (160, 64, "Deathmatch");

	M_PrintWhite (0, 72, "        Teamplay");
	if (rogue)
	{
		char *msg;

		switch((int)teamplay.value)
		{
			case 1: msg = "No Friendly Fire"; break;
			case 2: msg = "Friendly Fire"; break;
			case 3: msg = "Tag"; break;
			case 4: msg = "Capture the Flag"; break;
			case 5: msg = "One Flag CTF"; break;
			case 6: msg = "Three Team CTF"; break;
			default: msg = "Off"; break;
		}
		M_Print (160, 72, msg);
	}
	else
	{
		char *msg;

		switch((int)teamplay.value)
		{
			case 1: msg = "No Friendly Fire"; break;
			case 2: msg = "Friendly Fire"; break;
			default: msg = "Off"; break;
		}
		M_Print (160, 72, msg);
	}

	M_PrintWhite (0, 80, "            Skill");
	if (skill.value == 0)
		M_Print (160, 80, "Easy difficulty");
	else if (skill.value == 1)
		M_Print (160, 80, "Normal difficulty");
	else if (skill.value == 2)
		M_Print (160, 80, "Hard difficulty");
	else
    {
        if (kurok)
            M_Print (160, 80, "Insane difficulty");
        else
            M_Print (160, 80, "Nightmare difficulty");
    }

	M_PrintWhite (0, 88, "       Frag Limit");
	if (fraglimit.value == 0)
		M_Print (160, 88, "none");
	else
		M_Print (160, 88, va("%i frags", (int)fraglimit.value));

	M_PrintWhite (0, 96, "       Time Limit");
	if (timelimit.value == 0)
		M_Print (160, 96, "none");
	else
		M_Print (160, 96, va("%i minutes", (int)timelimit.value));

	M_PrintWhite (0, 104, "        Auto Aim");
	if (sv_aim.value == 1)
		M_Print (160, 104, "Off");
	else
		M_Print (160, 104, "On");

	M_PrintWhite (0, 112, "      Level Exits");
	if (noexit.value == 1)
		M_Print (160, 112, "Off");
	else
		M_Print (160, 112, "On");

	M_PrintWhite (0, 128, "         Episode");
   //MED 01/06/97 added hipnotic episodes
   if (hipnotic)
      M_Print (160, 128, hipnoticepisodes[startepisode].description);
   //PGM 01/07/97 added rogue episodes
   else if (rogue)
      M_Print (160, 128, rogueepisodes[startepisode].description);
   else if (kurok)
      M_Print (160, 128, kurokepisodes[startepisode].description);
   else
      M_Print (160, 128, episodes[startepisode].description);

	M_PrintWhite (0, 136, "           Level");
   //MED 01/06/97 added hipnotic episodes
   if (hipnotic)
   {
      M_Print (160, 136, hipnoticlevels[hipnoticepisodes[startepisode].firstLevel + startlevel].description);
      M_Print (160, 144, hipnoticlevels[hipnoticepisodes[startepisode].firstLevel + startlevel].name);
   }
   //PGM 01/07/97 added rogue episodes
   else if (rogue)
   {
      M_Print (160, 136, roguelevels[rogueepisodes[startepisode].firstLevel + startlevel].description);
      M_Print (160, 144, roguelevels[rogueepisodes[startepisode].firstLevel + startlevel].name);
   }
   else if (kurok)
   {
      M_Print (160, 136, kuroklevels[kurokepisodes[startepisode].firstLevel + startlevel].description);
      M_Print (160, 144, kuroklevels[kurokepisodes[startepisode].firstLevel + startlevel].name);
   }
   else
   {
      M_Print (160, 136, levels[episodes[startepisode].firstLevel + startlevel].description);
      M_Print (160, 144, levels[episodes[startepisode].firstLevel + startlevel].name);
   }
   
// line cursor
	M_DrawCharacter (144, gameoptions_cursor_table[gameoptions_cursor], 12+((int)(realtime*4)&1));

	if (m_serverInfoMessage)
	{
		if ((realtime - m_serverInfoMessageTime) < 5.0)
		{
			x = (320-26*8)/2;
			M_DrawTextBox (x, 138, 24, 4);
			x += 8;
			M_Print (x, 146, "   More than 4 players  ");
			M_Print (x, 154, " requires using command ");
			M_Print (x, 162, "   line -listen. Use    ");
			M_Print (x, 170, " -listen 8 for example. ");
		}
		else
		{
			m_serverInfoMessage = false;
		}
	}
}


void M_NetStart_Change (int dir)
{
	int count;

	switch (gameoptions_cursor)
	{
	case 1:
		maxplayers += dir;
		if (maxplayers > svs.maxclientslimit)
		{
			maxplayers = svs.maxclientslimit;
			m_serverInfoMessage = true;
			m_serverInfoMessageTime = realtime;
		}
		if (maxplayers < 2)
			maxplayers = 2;
		break;

	case 2:
		Cvar_SetValue ("coop", coop.value ? 0 : 1);
		break;

	case 3:
		if (rogue)
			count = 6;
		else
			count = 2;

		Cvar_SetValue ("teamplay", teamplay.value + dir);
		if (teamplay.value > count)
			Cvar_SetValue ("teamplay", 0);
		else if (teamplay.value < 0)
			Cvar_SetValue ("teamplay", count);
		break;

	case 4:
		Cvar_SetValue ("skill", skill.value + dir);
		if (skill.value > 3)
			Cvar_SetValue ("skill", 0);
		if (skill.value < 0)
			Cvar_SetValue ("skill", 3);
		break;

	case 5:
		Cvar_SetValue ("fraglimit", fraglimit.value + dir*10);
		if (fraglimit.value > 100)
			Cvar_SetValue ("fraglimit", 0);
		if (fraglimit.value < 0)
			Cvar_SetValue ("fraglimit", 100);
		break;

	case 6:
		Cvar_SetValue ("timelimit", timelimit.value + dir*5);
		if (timelimit.value > 60)
			Cvar_SetValue ("timelimit", 0);
		if (timelimit.value < 0)
			Cvar_SetValue ("timelimit", 60);
		break;

	case 7:
/*
		sv_aim.value += dir * 0.01;
		if (sv_aim.value < 0.9)
			sv_aim.value = 0.9;
		if (sv_aim.value > 1)
			sv_aim.value = 1;
*/
		Cvar_SetValue ("sv_aim", sv_aim.value + dir * 0.01);
		if (sv_aim.value > 1)
			Cvar_SetValue ("sv_aim", 0.99);
		if (sv_aim.value < 0.99)
			Cvar_SetValue ("sv_aim", 1);
		break;

	case 8:
		Cvar_SetValue ("noexit", noexit.value ? 0 : 1);
		break;

	case 9:
		startepisode += dir;
	//MED 01/06/97 added hipnotic count
		if (hipnotic)
			count = 6;
	//PGM 01/07/97 added rogue count
	//PGM 03/02/97 added 1 for dmatch episode
		else if (rogue)
			count = 4;
		else if (kurok)
			count = 3;
		else if (registered.value)
			count = 7;
		else
			count = 2;

		if (startepisode < 0)
			startepisode = count - 1;

		if (startepisode >= count)
			startepisode = 0;

		startlevel = 0;
		break;

	case 10:
		startlevel += dir;
    //MED 01/06/97 added hipnotic episodes
		if (hipnotic)
			count = hipnoticepisodes[startepisode].levels;
	//PGM 01/06/97 added hipnotic episodes
		else if (rogue)
			count = rogueepisodes[startepisode].levels;
		else if (kurok)
			count = kurokepisodes[startepisode].levels;
		else
			count = episodes[startepisode].levels;

		if (startlevel < 0)
			startlevel = count - 1;

		if (startlevel >= count)
			startlevel = 0;
		break;
	}
}

void M_GameOptions_Key (int key)
{
	switch (key)
	{
	case K_ESCAPE:
//		M_Menu_Net_f ();
		M_Menu_MultiPlayer_f ();
		break;

	case K_UPARROW:
		S_LocalSound ("misc/menu1.wav");
		gameoptions_cursor--;
		if (gameoptions_cursor < 0)
			gameoptions_cursor = NUM_GAMEOPTIONS-1;
		break;

	case K_DOWNARROW:
		S_LocalSound ("misc/menu1.wav");
		gameoptions_cursor++;
		if (gameoptions_cursor >= NUM_GAMEOPTIONS)
			gameoptions_cursor = 0;
		break;

	case K_LEFTARROW:
		if (gameoptions_cursor == 0)
			break;
		S_LocalSound ("misc/menu3.wav");
		M_NetStart_Change (-1);
		break;

	case K_RIGHTARROW:
		if (gameoptions_cursor == 0)
			break;
		S_LocalSound ("misc/menu3.wav");
		M_NetStart_Change (1);
		break;

	case K_ENTER:
		S_LocalSound ("misc/menu2.wav");
		if (gameoptions_cursor == 0)
		{
			if (sv.active)
				Cbuf_AddText ("disconnect\n");
			Cbuf_AddText ("listen 0\n");	// so host_netport will be re-examined
			Cbuf_AddText ( va ("maxplayers %u\n", maxplayers) );
			SCR_BeginLoadingPlaque ();

			if (hipnotic)
				Cbuf_AddText ( va ("map %s\n", hipnoticlevels[hipnoticepisodes[startepisode].firstLevel + startlevel].name) );
			else if (rogue)
				Cbuf_AddText ( va ("map %s\n", roguelevels[rogueepisodes[startepisode].firstLevel + startlevel].name) );
			else if (kurok)
				Cbuf_AddText ( va ("map %s\n", kuroklevels[kurokepisodes[startepisode].firstLevel + startlevel].name) );
			else
				Cbuf_AddText ( va ("map %s\n", levels[episodes[startepisode].firstLevel + startlevel].name) );

			return;
		}

		M_NetStart_Change (1);
		break;
	}
}

//=============================================================================
/* SEARCH MENU */

qboolean	searchComplete = false;
double		searchCompleteTime;

void M_Menu_Search_f (void)
{
	key_dest = key_menu;
	m_state = m_search;
	m_entersound = false;
	slistSilent = true;
	slistLocal = false;
	searchComplete = false;
	NET_Slist_f();

}


void M_Search_Draw (void)
{
	qpic_t	*p;
	int x;

	p = Draw_CachePic ("gfx/p_multi.lmp");
	M_DrawPic ( (320-p->width)/2, 4, p);
	x = (320/2) - ((12*8)/2) + 4;
	M_DrawTextBox (x-8, 32, 12, 1);
	M_Print (x, 40, "Searching...");

	if(slistInProgress)
	{
		NET_Poll();
		return;
	}

	if (! searchComplete)
	{
		searchComplete = true;
		searchCompleteTime = realtime;
	}

	if (hostCacheCount)
	{
		M_Menu_ServerList_f ();
		return;
	}

    if (kurok)
    {
		M_PrintWhite ((320/2) - ((22*8)/2), 64, "No Kurok servers found");
		if ((realtime - searchCompleteTime) < 3.0)
			return;
    }
    
    else
    {
		M_PrintWhite ((320/2) - ((22*8)/2), 64, "No servers found");
		if ((realtime - searchCompleteTime) < 3.0)
			return;
    }
    
	M_Menu_LanConfig_f ();
}


void M_Search_Key (int key)
{
}

//=============================================================================
/* SLIST MENU */

int		slist_cursor;
qboolean slist_sorted;

void M_Menu_ServerList_f (void)
{
	key_dest = key_menu;
	m_state = m_slist;
	m_entersound = true;
	slist_cursor = 0;
	m_return_onerror = false;
	m_return_reason[0] = 0;
	slist_sorted = false;
}


void M_ServerList_Draw (void)
{
	int		n;
	char	string [64];
	qpic_t	*p;

	if (!slist_sorted)
	{
		if (hostCacheCount > 1)
		{
			int	i,j;
			hostcache_t temp;
			for (i = 0; i < hostCacheCount; i++)
				for (j = i+1; j < hostCacheCount; j++)
					if (strcmp(hostcache[j].name, hostcache[i].name) < 0)
					{
						Q_memcpy(&temp, &hostcache[j], sizeof(hostcache_t));
						Q_memcpy(&hostcache[j], &hostcache[i], sizeof(hostcache_t));
						Q_memcpy(&hostcache[i], &temp, sizeof(hostcache_t));
					}
		}
		slist_sorted = true;
	}

    if(kurok)
        M_DrawCharacter (0, 32 + slist_cursor*8, 12+((int)(realtime*30)&1));
    else
        M_DrawCharacter (0, 32 + slist_cursor*8, 12+((int)(realtime*4)&1));

	p = Draw_CachePic ("gfx/p_multi.lmp");
	M_DrawPic ( (320-p->width)/2, 4, p);
	for (n = 0; n < hostCacheCount; n++)
	{
		if (hostcache[n].maxusers)
			sprintf(string, "%-15.15s %-15.15s %2u/%2u\n", hostcache[n].name, hostcache[n].map, hostcache[n].users, hostcache[n].maxusers);
		else
			sprintf(string, "%-15.15s %-15.15s\n", hostcache[n].name, hostcache[n].map);
		M_Print (16, 32 + 8*n, string);
	}

	if (*m_return_reason)
		M_PrintWhite (16, 148, m_return_reason);
}


void M_ServerList_Key (int k)
{
	switch (k)
	{
	case K_ESCAPE:
		M_Menu_LanConfig_f ();
		break;

	case K_SPACE:
		M_Menu_Search_f ();
		break;

	case K_UPARROW:
	case K_LEFTARROW:
		S_LocalSound ("misc/menu1.wav");
		slist_cursor--;
		if (slist_cursor < 0)
			slist_cursor = hostCacheCount - 1;
		break;

	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound ("misc/menu1.wav");
		slist_cursor++;
		if (slist_cursor >= hostCacheCount)
			slist_cursor = 0;
		break;

	case K_ENTER:
		S_LocalSound ("misc/menu2.wav");
		m_return_state = m_state;
		m_return_onerror = true;
		slist_sorted = false;
		key_dest = key_game;
		m_state = m_none;
		Cbuf_AddText ( va ("connect \"%s\"\n", hostcache[slist_cursor].cname) );
		break;

	default:
		break;
	}

}

//=============================================================================
/* Menu Subsystem */


void M_Init (void)
{
	Cmd_AddCommand ("togglemenu", M_ToggleMenu_f);

	Cmd_AddCommand ("menu_main", M_Menu_Main_f);
	Cmd_AddCommand ("menu_singleplayer", M_Menu_SinglePlayer_f);
	Cmd_AddCommand ("menu_load", M_Menu_Load_f);
	Cmd_AddCommand ("menu_save", M_Menu_Save_f);
	Cmd_AddCommand ("menu_multiplayer", M_Menu_MultiPlayer_f);
	Cmd_AddCommand ("menu_setup", M_Menu_Setup_f);
	Cmd_AddCommand ("menu_options", M_Menu_Options_f);
	Cmd_AddCommand ("menu_keys", M_Menu_Keys_f);
	Cmd_AddCommand ("menu_video", M_Menu_Video_f);
	Cmd_AddCommand ("help", M_Menu_Help_f);
	Cmd_AddCommand ("menu_quit", M_Menu_Quit_f);
}


void M_Draw (void)
{
	if (m_state == m_none || key_dest != key_menu)
		return;

	if (!m_recursiveDraw)
	{
		scr_copyeverything = 1;

		if (scr_con_current)
		{
			Draw_ConsoleBackground (vid.height);
			VID_UnlockBuffer ();
			S_ExtraUpdate ();
			VID_LockBuffer ();
		}
		else
		{
            if (kurok)
                Draw_FadeScreen2 ();
            else
			    Draw_FadeScreen ();
        }
		scr_fullupdate = 0;
	}
	else
	{
		m_recursiveDraw = false;
	}

	switch (m_state)
	{
	case m_none:
		break;

	case m_main:
		M_Main_Draw ();
		break;

	case m_singleplayer:
		M_SinglePlayer_Draw ();
		break;

	case m_load:
		M_Load_Draw ();
		break;

	case m_save:
		M_Save_Draw ();
		break;

	case m_multiplayer:
		M_MultiPlayer_Draw ();
		break;

	case m_setup:
		M_Setup_Draw ();
		break;

	case m_net:
		M_Net_Draw ();
		break;

	case m_options:
		M_Options_Draw ();
		break;

	case m_keys:
		M_Keys_Draw ();
		break;

	case m_video:
		M_Video_Draw ();
		break;

	case m_help:
		M_Help_Draw ();
		break;

	case m_quit:
		M_Quit_Draw ();
		break;

	case m_serialconfig:
		M_SerialConfig_Draw ();
		break;

	case m_modemconfig:
		M_ModemConfig_Draw ();
		break;

	case m_lanconfig:
		M_LanConfig_Draw ();
		break;

	case m_gameoptions:
		M_GameOptions_Draw ();
		break;

	case m_search:
		M_Search_Draw ();
		break;

	case m_slist:
		M_ServerList_Draw ();
		break;
		
	case m_osk:
		M_OSK_Draw();
		break;
	}

	if (m_entersound)
	{
		S_LocalSound ("misc/menu2.wav");
		m_entersound = false;
	}

	VID_UnlockBuffer ();
	S_ExtraUpdate ();
	VID_LockBuffer ();
}


void M_Keydown (int key)
{
	switch (m_state)
	{
	case m_none:
		return;

	case m_main:
		M_Main_Key (key);
		return;

	case m_singleplayer:
		M_SinglePlayer_Key (key);
		return;

	case m_load:
		M_Load_Key (key);
		return;

	case m_save:
		M_Save_Key (key);
		return;

	case m_multiplayer:
		M_MultiPlayer_Key (key);
		return;

	case m_setup:
		M_Setup_Key (key);
		return;

	case m_net:
		M_Net_Key (key);
		return;

	case m_options:
		M_Options_Key (key);
		return;

	case m_keys:
		M_Keys_Key (key);
		return;

	case m_video:
		M_Video_Key (key);
		return;

	case m_help:
		M_Help_Key (key);
		return;

	case m_quit:
		M_Quit_Key (key);
		return;

	case m_serialconfig:
		M_SerialConfig_Key (key);
		return;

	case m_modemconfig:
		M_ModemConfig_Key (key);
		return;

	case m_lanconfig:
		M_LanConfig_Key (key);
		return;

	case m_gameoptions:
		M_GameOptions_Key (key);
		return;

	case m_search:
		M_Search_Key (key);
		break;
			
	case m_slist:
		M_ServerList_Key (key);
		return;

	case m_osk:
		M_OSK_Key(key);	
	}
}


void M_ConfigureNetSubsystem(void)
{
// enable/disable net systems to match desired config

	Cbuf_AddText ("stopdemo\n");
	if (SerialConfig || DirectConfig)
	{
		Cbuf_AddText ("com1 enable\n");
	}

	if (IPXConfig || TCPIPConfig)
		net_hostport = lanConfig_port;
}
