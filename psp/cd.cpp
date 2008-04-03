/*
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

#include <cstddef>
#include <stdio.h>

#include <pspaudiolib.h>

//#include <mad.h>

#include "mp3.h"

extern "C"
{
#include "../quakedef.h"
}

extern 	cvar_t bgmtype;
extern	cvar_t bgmvolume;


namespace quake
{
	namespace cd
	{
		struct Sample
		{
			short left;
			short right;
		};

//		static int				file = -1;
		static int		        last_track = 0;

		static bool	 playing  = false;
		static bool	 paused   = false;
		static bool	 enabled  = false;
		static float cdvolume = 0;
	}
}

using namespace quake;
using namespace quake::cd;

static void CD_f (void)
{
	char	*command;

	if (Cmd_Argc() < 2)
	{
		Con_Printf("commands:");
		Con_Printf("on, off, reset, remap, \n");
		Con_Printf("play, stop, loop, pause, resume\n");
		Con_Printf("eject, close, info\n");
		return;
	}

	command = Cmd_Argv (1);

	if (Q_strcasecmp(command, "on") == 0)
	{
		enabled = true;
		return;
	}

	if (Q_strcasecmp(command, "off") == 0)
	{
		if (playing)
			CDAudio_Stop();
		enabled = false;
		return;
	}

	if (Q_strcasecmp(command, "reset") == 0)
	{
		enabled = true;
		if (playing)
			CDAudio_Stop();
		return;
	}

	if (Q_strcasecmp(command, "remap") == 0)
	{
		return;
	}

	if (Q_strcasecmp(command, "close") == 0)
	{
		return;
	}

	if (Q_strcasecmp(command, "play") == 0)
	{
		CDAudio_Play((byte)atoi(Cmd_Argv (2)), (qboolean) false);
		return;
	}

	if (Q_strcasecmp(command, "loop") == 0)
	{
		CDAudio_Play((byte)atoi(Cmd_Argv (2)), (qboolean) true);
		return;
	}

	if (Q_strcasecmp(command, "stop") == 0)
	{
		CDAudio_Stop();
		return;
	}

	if (Q_strcasecmp(command, "pause") == 0)
	{
		CDAudio_Pause();
		return;
	}

	if (Q_strcasecmp(command, "resume") == 0)
	{
		CDAudio_Resume();
		return;
	}

	if (Q_strcasecmp(command, "eject") == 0)
	{
		if (playing)
			CDAudio_Stop();
		return;
	}

	if (Q_strcasecmp(command, "info") == 0)
	{
		return;
	}
}

void CDAudio_VolumeChange(float bgmvolume)
{
	int volume = (int) (bgmvolume * (float) PSP_VOLUME_MAX);
	//pspAudioSetVolume(1, volume, volume);
	mp3_volume = volume;
	cdvolume = bgmvolume;
	
	changeMp3Volume=0;
	
//	Con_Printf("Volume changed to : %i\n", mp3_volume);
	
}

extern "C" int sceKernelDelayThread(int delay); 

void CDAudio_Play(byte track, qboolean looping)
{
	last_track = track;
	CDAudio_Stop();

	char path[256];
	sprintf(path, "%s\\MP3\\%02u.mp3", host_parms.basedir, track);
	
	int ret;
	ret = mp3_start_play(path, 0);

	if(ret != 2)
	{
		Con_Printf("Playing %s\n", path);
		playing = true;
	}
	else
	{
		Con_Printf("Couldn't find %s\n", path);
		playing = false;
		CDAudio_VolumeChange(0);
	}
	

	CDAudio_VolumeChange(bgmvolume.value);
}

void CDAudio_Stop(void)
{
	mp3_job_started = 0;
		
//	file = -1;
	playing = false;
	CDAudio_VolumeChange(0);
}

void CDAudio_Pause(void)
{
	
	CDAudio_VolumeChange(0);
	paused = true;
}

void CDAudio_Resume(void)
{
	
	CDAudio_VolumeChange(bgmvolume.value);
	paused = false;
}

void CDAudio_Update(void)
{
	
	//if (bgmvolume.value != mp3_volume)
	//	CDAudio_VolumeChange(bgmvolume.value);
	//if(changeMp3Volume) CDAudio_VolumeChange(bgmvolume.value);

	if (strcmpi(bgmtype.string,"cd") == 0) {
		if (playing == false) {
			CDAudio_Play(last_track, (qboolean) false);
		}
		if (paused == true) {
			CDAudio_Resume();
		}
	
	} else {
		if (paused == false) {
			CDAudio_Pause();
		}
		if (playing == true) {
			CDAudio_Stop();
		}
	}
}

int CDAudio_Init(void)
{
		
	mp3_init();
	sceKernelDelayThread(5*10000);
	
	enabled = true;
	
	Cmd_AddCommand ("cd", CD_f);
	
	Con_Printf("CD Audio Initialized\n");

	return 0;
}

void CDAudio_Shutdown(void)
{
	Con_Printf("CDAudio_Shutdown\n");
	
	CDAudio_Stop();

	sceKernelDelayThread(5*10000);
//	mp3_deinit();

}

