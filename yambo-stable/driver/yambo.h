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
 int yambo_driver
#else
 int yambo_driver_
#endif
 (int *, int *,int *,int *,int *,int *,int *,int *,
  char *rnstr2, char *inf, char *id, char *od, char *com_dir, char *js,
  int lni,int iif,int iid,int iod,int icd,int ijs);
/*
 Command line structure
*/
 static Ldes opts[] = { /* Int Real Ch (dummy) Parallel_option*/
  {"help",   "h","Short Help",0,0,0,0,0},              
  {"lhelp",  "H","Long Help",0,0,0,0,0},               
  {"jobstr", "J","Job string identifier",0,0,1,0,1},   
  {"infver", "V","Input file verbosity",0,0,1,0,0},    
  {"DESC",   " ","[opt=RL,kpt,sc,qp,io,gen,resp,all,par]",0,0,0,0,0},
  {"ifile",  "F","Input file",0,0,1,0,1},              
  {"idir",   "I","Core I/O directory",0,0,1,0,1},         
  {"odir",   "O","Additional I/O directory",0,0,1,0,1},   
  {"cdir",   "C","Communications I/O directory",0,0,1,0,1},
  {"dbpr",   "D","DataBases properties",0,0,0,0,0},    
  {"wallt",  "W","Wall Time limitation (1d2h30m format)",0,0,1,0,1}, 
  {"quiet",  "Q","Don't launch the text editor",0,0,0,0,0}, 
#if defined _MPI
  {"parenv", "E","Environment Parallel Variables file",0,0,1,0,1},               
  {"nompi",  "M","Switch-off MPI support (serial run)",0,0,0,0,0}, 
#endif
#if defined _OPENMP
  {"noopenmp","N","Switch-off OpenMP support (single thread run)",0,0,0,0,1}, 
#endif
  {"setup",  "i","Initialization",0,0,0,0,0},          
  {"optics", "o","Optics [opt=(c)hi is (G)-space / (b)se is (eh)-space ]",0,0,1,0,0},
  {"kernel", "k","Kernel [opt=hartree/alda/lrc/hf/sex]",0,0,1,0,0},
  {"DESC",   " ","       (hf/sex only eh-space; lrc only G-space)",0,0,0,0,0},
  {"bss",    "y","BSE solver [opt=h/d/(p/f)i]",0,0,1,0,0},                      
  {"DESC",   " ","       (h)aydock/(d)iagonalization/(i)nversion",0,0,0,0,0},
  {"rim_cut","r","Coulomb potential",0,0,0,0,0},  
  {"HF_and_locXC",  "x","Hartree-Fock Self-energy and local XC",0,0,0,0,0},      
  {"em1d",   "d","Dynamical Inverse Dielectric Matrix",0,0,0,0,0},     
  {"em1s",   "b","Static Inverse Dielectric Matrix",0,0,0,0,0},        
  {"gwapprx","p","GW approximations [opt=(p)PA/(c)HOSEX]",0,0,1,0,0},              
  {"gw0",    "g","Dyson Equation solver",0,0,1,0,0}, 
  {"DESC",   " ","[opt=(n)ewton/(s)ecant/(g)reen]",0,0,0,0,0},
  {"life",   "l","GoWo Quasiparticle lifetimes",0,0,0,0,0},                  
  {"acfdt",  "a","ACFDT Total Energy",0,0,0,0,0},                            
#if defined _RT
  {"tdBSE", "q","Real-time Bethe-Salpter Equation ",0,0,0,0,0}, 
#endif
#if defined _ELPH 
  {"corrtp", "c","Correlation [opt=e(lectrons)/p(honons)/b(oth)]",0,0,1,0,0},    
  {"ElPhHam","f","Electron-Phonon Hamiltonian",0,0,0,0,0},    
#endif
#if defined _RT
  {"collisions", "e","Evaluate Collisions",0,0,0,0,0}, 
#endif
  {NULL,NULL,NULL,0,0,0,0,0}
 };
 char *tool="yambo";
 char *tdesc="A shiny pot of fun and happiness [C.D.Hogan]";
