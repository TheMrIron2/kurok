#include "quakedef.h"
#include "alpha.h"

float        model_alpha;
int 		 eval_alpha;

int FindFieldOffset(char *field)
{
	ddef_t *d;
	d = ED_FindField(field);
	if (!d)
		return 0;
	return d->ofs*4;
}

void FindEdictFieldOffsets() {
//        dfunction_t *f;

//	eval_gravity = FindFieldOffset("gravity");
	eval_alpha = FindFieldOffset("alpha");
//	eval_fullbright = FindFieldOffset("fullbright");
//	eval_idealpitch = FindFieldOffset("idealpitch");
//	eval_pitch_speed = FindFieldOffset("pitch_speed");

//        RestoreGame = 0;
//        if ((f = ED_FindFunction ("RestoreGame")) != NULL)
//           RestoreGame = (func_t)(f - pr_functions);
};
