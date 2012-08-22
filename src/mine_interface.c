/* This code is written by Michele Filosi  <michele.filosi@gmail.com>
   Roberto Visintainer <r.visintainer@gmail.com>.
   2012 

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdlib.h>
#include <math.h>
#include "core.h"
#include "mine.h"
#include <R.h>
#include <Rinternals.h>

SEXP mineRonevar (SEXP x, SEXP y, SEXP alpha, SEXP C){
  
  double *restmp;
  mine_problem *prob;
  mine_parameter *param;
  mine_score *minescore;
  SEXP res;
  
  PROTECT(alpha = coerceVector(alpha,REALSXP));
  PROTECT(C = coerceVector(C,INTSXP));
  PROTECT(res=allocVector(REALSXP,4));
  restmp=REAL(res);
    
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=asReal(alpha);
  param->c=asReal(C);
  
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=length(x);
  prob->x=REAL(x);
  prob->y=REAL(y);
  
  minescore=mine_compute_score(prob,param);
  restmp[0]=mic(minescore);
  restmp[1]=mas(minescore);
  restmp[2]=mev(minescore);
  restmp[3]=mcn(minescore);
  
  /* Free */
  Free(prob);
  Free(param);
  mine_free_score(&minescore);
  UNPROTECT(3);
  return(res);
}




SEXP mineRall (SEXP x, SEXP nrx, SEXP ncx, SEXP alpha, SEXP C)
{
  R_len_t i, j, rx, cx;
  double score;
  double **pointers;
  mine_problem *prob;
  mine_parameter *param;
  mine_score *minescore;
  SEXP res, mydim, resmic, resmas, resmev, resmcn, names;
  
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=asReal(alpha);
  param->c=asReal(C);
  
  /* Matrix dimension */
  rx=asInteger(nrx);
  cx=asInteger(ncx);
  
  PROTECT(x=coerceVector(x,REALSXP));
  
  /* Initialize data matrix */
  pointers = (double **) Calloc(cx, double *);
  for (i = 0; i<cx;i++){
    pointers[i] = (double *) Calloc(rx, double);
    pointers[i] = &REAL(x)[i * rx];
  }
  
  /* Initialize result matrix */
  PROTECT(resmic=allocVector(REALSXP,cx*cx));
  PROTECT(resmas=allocVector(REALSXP,cx*cx));
  PROTECT(resmev=allocVector(REALSXP,cx*cx));
  PROTECT(resmcn=allocVector(REALSXP,cx*cx));
  PROTECT(res=allocVector(VECSXP,4));
  
  /* Allocating result list */
  SET_VECTOR_ELT(res, 0, resmic);
  SET_VECTOR_ELT(res, 1, resmas);
  SET_VECTOR_ELT(res, 2, resmev);
  SET_VECTOR_ELT(res, 3, resmcn);
  
  /* Set the mine_problem */
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=rx;
  
  for (i = 0; i<cx; i++){
    prob->x=pointers[i];
    for (j = i; j<cx; j++){
      prob->y=pointers[j];
      /* Computing MINE scores */
      minescore=mine_compute_score(prob,param);
      score=mic(minescore);
      REAL(resmic)[(cx*j) + i] = score;
      REAL(resmic)[(cx*i) + j] = score;

      score=mas(minescore);
      REAL(resmas)[(cx*j) + i] = score;
      REAL(resmas)[(cx*i) + j] = score;

      score=mev(minescore);
      REAL(resmev)[(cx*j) + i] = score;
      REAL(resmev)[(cx*i) + j] = score;

      score=mcn(minescore);
      REAL(resmcn)[(cx*j) + i] = score;
      REAL(resmcn)[(cx*i) + j] = score;
      
      /* Free score */
      mine_free_score(&minescore);
    }
  }
    
  /* Set matrix dimension to be passed to R object */
  PROTECT(mydim=allocVector(INTSXP, 2));
  INTEGER(mydim)[0] = cx;
  INTEGER(mydim)[1] = cx;
  
  setAttrib(resmic, R_DimSymbol, mydim);
  setAttrib(resmas, R_DimSymbol, mydim);
  setAttrib(resmev, R_DimSymbol, mydim);
  setAttrib(resmcn, R_DimSymbol, mydim);

  PROTECT(names=allocVector(STRSXP,4));
  SET_STRING_ELT(names,0,mkChar("MIC"));
  SET_STRING_ELT(names,1,mkChar("MAS"));  
  SET_STRING_ELT(names,2,mkChar("MEV"));  
  SET_STRING_ELT(names,3,mkChar("MCN"));
  setAttrib(res, R_NamesSymbol,names);
  
  /* Free memeory */
  UNPROTECT(8);
  Free(pointers);
  Free(param);
  Free(prob);
  /* Return a named list of 4 matrix */
  return(res);
}

