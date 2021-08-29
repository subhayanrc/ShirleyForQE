/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Mon Jan 16 09:08:03 EST 2017 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 5 -name t2_5 -include t.h */

/*
 * This function contains 44 FP additions, 40 FP multiplications,
 * (or, 14 additions, 10 multiplications, 30 fused multiply/add),
 * 47 stack variables, 4 constants, and 20 memory accesses
 */
#include "t.h"

static void t2_5(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP618033988, +0.618033988749894848204586834365638117720309180);
     {
	  INT m;
	  for (m = mb, W = W + (mb * 4); m < me; m = m + 1, ri = ri + ms, ii = ii + ms, W = W + 4, MAKE_VOLATILE_STRIDE(10, rs)) {
	       E Ta, T1, TO, Tp, TS, Ti, TL, TC, To, TE, Ts, TF, T2, T8, T5;
	       E TT, Tt, TG;
	       T2 = W[0];
	       Ta = W[3];
	       T8 = W[2];
	       T5 = W[1];
	       {
		    E Tq, Tr, Te, T9;
		    T1 = ri[0];
		    Te = T2 * Ta;
		    T9 = T2 * T8;
		    TO = ii[0];
		    {
			 E T3, Tf, Tm, Tj, Tb, T4, T6, Tc, Tg;
			 T3 = ri[WS(rs, 1)];
			 Tf = FMA(T5, T8, Te);
			 Tm = FNMS(T5, T8, Te);
			 Tj = FMA(T5, Ta, T9);
			 Tb = FNMS(T5, Ta, T9);
			 T4 = T2 * T3;
			 T6 = ii[WS(rs, 1)];
			 Tc = ri[WS(rs, 4)];
			 Tg = ii[WS(rs, 4)];
			 {
			      E Tk, Tl, Tn, TD;
			      {
				   E T7, Tz, Th, TB, Ty, Td, TA;
				   Tk = ri[WS(rs, 2)];
				   T7 = FMA(T5, T6, T4);
				   Ty = T2 * T6;
				   Td = Tb * Tc;
				   TA = Tb * Tg;
				   Tl = Tj * Tk;
				   Tz = FNMS(T5, T3, Ty);
				   Th = FMA(Tf, Tg, Td);
				   TB = FNMS(Tf, Tc, TA);
				   Tn = ii[WS(rs, 2)];
				   Tp = ri[WS(rs, 3)];
				   TS = T7 - Th;
				   Ti = T7 + Th;
				   TL = Tz + TB;
				   TC = Tz - TB;
				   TD = Tj * Tn;
				   Tq = T8 * Tp;
				   Tr = ii[WS(rs, 3)];
			      }
			      To = FMA(Tm, Tn, Tl);
			      TE = FNMS(Tm, Tk, TD);
			 }
		    }
		    Ts = FMA(Ta, Tr, Tq);
		    TF = T8 * Tr;
	       }
	       TT = To - Ts;
	       Tt = To + Ts;
	       TG = FNMS(Ta, Tp, TF);
	       {
		    E TU, TW, TV, TR, Tw, Tu;
		    TU = FMA(KP618033988, TT, TS);
		    TW = FNMS(KP618033988, TS, TT);
		    Tw = Ti - Tt;
		    Tu = Ti + Tt;
		    {
			 E TM, TH, Tv, TI, TK;
			 TM = TE + TG;
			 TH = TE - TG;
			 ri[0] = T1 + Tu;
			 Tv = FNMS(KP250000000, Tu, T1);
			 TI = FMA(KP618033988, TH, TC);
			 TK = FNMS(KP618033988, TC, TH);
			 {
			      E TQ, TN, TJ, Tx, TP;
			      TQ = TL - TM;
			      TN = TL + TM;
			      TJ = FNMS(KP559016994, Tw, Tv);
			      Tx = FMA(KP559016994, Tw, Tv);
			      ii[0] = TN + TO;
			      TP = FNMS(KP250000000, TN, TO);
			      ri[WS(rs, 1)] = FMA(KP951056516, TI, Tx);
			      ri[WS(rs, 4)] = FNMS(KP951056516, TI, Tx);
			      ri[WS(rs, 3)] = FMA(KP951056516, TK, TJ);
			      ri[WS(rs, 2)] = FNMS(KP951056516, TK, TJ);
			      TV = FNMS(KP559016994, TQ, TP);
			      TR = FMA(KP559016994, TQ, TP);
			 }
		    }
		    ii[WS(rs, 4)] = FMA(KP951056516, TU, TR);
		    ii[WS(rs, 1)] = FNMS(KP951056516, TU, TR);
		    ii[WS(rs, 3)] = FNMS(KP951056516, TW, TV);
		    ii[WS(rs, 2)] = FMA(KP951056516, TW, TV);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 0, 1},
     {TW_CEXP, 0, 3},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 5, "t2_5", twinstr, &GENUS, {14, 10, 30, 0}, 0, 0, 0 };

void X(codelet_t2_5) (planner *p) {
     X(kdft_dit_register) (p, t2_5, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle.native -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 5 -name t2_5 -include t.h */

/*
 * This function contains 44 FP additions, 32 FP multiplications,
 * (or, 30 additions, 18 multiplications, 14 fused multiply/add),
 * 37 stack variables, 4 constants, and 20 memory accesses
 */
#include "t.h"

static void t2_5(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     {
	  INT m;
	  for (m = mb, W = W + (mb * 4); m < me; m = m + 1, ri = ri + ms, ii = ii + ms, W = W + 4, MAKE_VOLATILE_STRIDE(10, rs)) {
	       E T2, T4, T7, T9, Tb, Tl, Tf, Tj;
	       {
		    E T8, Te, Ta, Td;
		    T2 = W[0];
		    T4 = W[1];
		    T7 = W[2];
		    T9 = W[3];
		    T8 = T2 * T7;
		    Te = T4 * T7;
		    Ta = T4 * T9;
		    Td = T2 * T9;
		    Tb = T8 - Ta;
		    Tl = Td - Te;
		    Tf = Td + Te;
		    Tj = T8 + Ta;
	       }
	       {
		    E T1, TI, Ty, TB, TN, TM, TF, TG, TH, Ti, Tr, Ts;
		    T1 = ri[0];
		    TI = ii[0];
		    {
			 E T6, Tw, Tq, TA, Th, Tx, Tn, Tz;
			 {
			      E T3, T5, To, Tp;
			      T3 = ri[WS(rs, 1)];
			      T5 = ii[WS(rs, 1)];
			      T6 = FMA(T2, T3, T4 * T5);
			      Tw = FNMS(T4, T3, T2 * T5);
			      To = ri[WS(rs, 3)];
			      Tp = ii[WS(rs, 3)];
			      Tq = FMA(T7, To, T9 * Tp);
			      TA = FNMS(T9, To, T7 * Tp);
			 }
			 {
			      E Tc, Tg, Tk, Tm;
			      Tc = ri[WS(rs, 4)];
			      Tg = ii[WS(rs, 4)];
			      Th = FMA(Tb, Tc, Tf * Tg);
			      Tx = FNMS(Tf, Tc, Tb * Tg);
			      Tk = ri[WS(rs, 2)];
			      Tm = ii[WS(rs, 2)];
			      Tn = FMA(Tj, Tk, Tl * Tm);
			      Tz = FNMS(Tl, Tk, Tj * Tm);
			 }
			 Ty = Tw - Tx;
			 TB = Tz - TA;
			 TN = Tn - Tq;
			 TM = T6 - Th;
			 TF = Tw + Tx;
			 TG = Tz + TA;
			 TH = TF + TG;
			 Ti = T6 + Th;
			 Tr = Tn + Tq;
			 Ts = Ti + Tr;
		    }
		    ri[0] = T1 + Ts;
		    ii[0] = TH + TI;
		    {
			 E TC, TE, Tv, TD, Tt, Tu;
			 TC = FMA(KP951056516, Ty, KP587785252 * TB);
			 TE = FNMS(KP587785252, Ty, KP951056516 * TB);
			 Tt = KP559016994 * (Ti - Tr);
			 Tu = FNMS(KP250000000, Ts, T1);
			 Tv = Tt + Tu;
			 TD = Tu - Tt;
			 ri[WS(rs, 4)] = Tv - TC;
			 ri[WS(rs, 3)] = TD + TE;
			 ri[WS(rs, 1)] = Tv + TC;
			 ri[WS(rs, 2)] = TD - TE;
		    }
		    {
			 E TO, TP, TL, TQ, TJ, TK;
			 TO = FMA(KP951056516, TM, KP587785252 * TN);
			 TP = FNMS(KP587785252, TM, KP951056516 * TN);
			 TJ = KP559016994 * (TF - TG);
			 TK = FNMS(KP250000000, TH, TI);
			 TL = TJ + TK;
			 TQ = TK - TJ;
			 ii[WS(rs, 1)] = TL - TO;
			 ii[WS(rs, 3)] = TQ - TP;
			 ii[WS(rs, 4)] = TO + TL;
			 ii[WS(rs, 2)] = TP + TQ;
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 0, 1},
     {TW_CEXP, 0, 3},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 5, "t2_5", twinstr, &GENUS, {30, 18, 14, 0}, 0, 0, 0 };

void X(codelet_t2_5) (planner *p) {
     X(kdft_dit_register) (p, t2_5, &desc);
}
#endif				/* HAVE_FMA */
