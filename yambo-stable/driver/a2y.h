/*
         Copyright (C) 2000-2017 the YAMBO team
               http://www.yambo-code.org
 
  Authors (see AUTHORS file for details): AM
  
  This file is distributed under the terms of the GNU 
  General Public License. You can redistribute it and/or 
  modify it under the terms of the GNU General Public 
  License as published by the Free Software Foundation; 
  either version 2, or (at your option) any later version.
 
  This program is distributed in the hope that it will 
  be useful, but WITHOUT ANY WARRANTY; without even the 
  implied warranty of MERCHANTABILITY or FITNESS FOR A 
  PARTICULAR PURPOSE.  See the GNU General Public License 
  for more details.
 
  You should have received a copy of the GNU General Public 
  License along with this program; if not, write to the Free 
  Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
  MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
*/
/*
 Driver declaration
*/
#if defined _FORTRAN_US
 int a2y_i_
#else
 int a2y_i
#endif
 (int *, int *,int *,int *,int *,int *,int *,int *,
  char *rnstr2, char *inf, char *id, char *od, char *com_dir, char *js,
  int lni,int iif,int iid,int iod,int icd,int ijs);
/*
 Command line structure
*/
 static Ldes opts[] = { /* Int Real Ch (dummy) Parallel_option*/
  {"help",  "h","Short Help",0,0,0,0,0}, 
  {"lhelp", "H","Long Help",0,0,0,0,0}, 
  {"nompi", "N","Skip MPI initialization",0,0,0,0,0},
  {"ifile", "F","Abinit file name/prefix",0,0,1,0,0}, 
  {"odir",  "O","Output directory",0,0,1,0,0}, 
  {"nodbfr","U","Do not fragment the DataBases (only for serial runs)",0,0,0,0,0},
  {"fragnb","b","Number of bands for each fragment",1,0,0,0,1},
  {"alat_f","a","Lattice constants rescaling factor",0,1,0,0,0},
  {"dupl",  "d","States duplication (artificial spin polarization)",0,0,0,0,0},
  {"notr",  "t","Force use of spatial Inv. instead of Time Rev.",0,0,0,0,0},
  {"nosy",  "n","Force no symmetries",0,0,0,0,0},
  {"nowf",  "w","Force no wavefunctions",0,0,0,0,0},
  {NULL,NULL,NULL,0,0,0,0,0}
 };
 char *tool="a2y";
 char *tdesc="A(binit) 2 Y(ambo) interface";
