/*  
    This code is written by Davide Albanese <davide.albanese@gmail.com>.
    (C) 2012 Davide Albanese, (C) 2012 Fondazione Bruno Kessler.

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
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "core.h"

#define MAX(a, b) ((a) > (b) ? (a):(b))
#define MIN(a, b) ((a) < (b) ? (a):(b))


/* HQ
 * Returns the entropy induced by the points on the 
 * partition Q.
 * See section 3.2.1 at page 10 in SOM.
 *
 * Parameters
 *   N : q x p matrix containing the number of points
 *       for each cell of the grid formed by P and Q
 *       partitions.
 *   q : number of rows of N (number of partitions in Q)
 *   p : number of cols of N (number of partitions in P)
 *   n : total number of points
 */
double HQ(int **N, int q, int p, int n)
{
  int i, j, sum;
  double prob, logprob, H;

  
  H = 0.0;
  for (i=0; i<q; i++)
    {
      sum = 0;
      for (j=0; j<p; j++)
	sum += N[i][j];

      prob = (double) sum / (double) n;
      if (prob != 0)
	{
	  logprob = log(prob);
	  H += prob * logprob;
	}
    }

  return -H;
}


/* HP3
 * Returns the entropy induced by the points on the 
 * partition <c_0, c_s, c_t>.
 * See line 5 of Algorithm 2 in SOM.
 *
 * Parameters
 *   np : p vector containing the number of points
 *        for each cell of the grid formed by P
 *   p : length of np (number of partitions in P)
 *   s : s in c_s
 *   t : t in c_t
 */
double HP3(int *np, int p, int s, int t)
{
  int i;
  int sum1, sum2, tot;
  double prob1, prob2, H;
  

  if (s == t)
    return 0.0;

  sum1 = 0;
  for (i=0; i<s; i++)
    sum1 += np[i];
  
  sum2 = 0;
  for (i=s; i<t; i++)
    sum2 += np[i];
  
  tot = sum1 + sum2;
  prob1 = (double) sum1 / (double) tot;
  prob2 = (double) sum2 / (double) tot;
  
  H = 0.0;
  
  if (prob1 != 0)
    H += prob1 * log(prob1);
  
  if (prob2 != 0)
    H += prob2 * log(prob2);

  return -H;
}


/* HPQ3
 * Returns the entropy induced by the points on the 
 * partition <c_0, c_s, c_t> and Q.
 * See line 5 of Algorithm 2 in SOM.
 *
 * Parameters
 *   N : q x p matrix containing the number of points
 *       for each cell of the grid formed by P and Q
 *       partitions.
 *   np : p vector containing the number of points
 *        for each cell of the grid formed by P
 *   q : number of rows of N (number of partitions in Q)
 *   p : number of cols of N and length of np (number of partitions in P)
 *   s : s in c_s
 *   t : t in c_t
 */
double HPQ3(int **N, int *np, int q, int p, int s, int t)
{
  int i, j;
  int sum1, sum2, tot;
  double prob1, prob2, H;
  

  tot = 0;
  for (i=0; i<t; i++)
    tot += np[i];

  H = 0.0;
  for (i=0; i<q; i++)
    {
      sum1 = 0;
      for (j=0; j<s; j++)
	sum1 += N[i][j];
      
      sum2 = 0;
      for (j=s; j<t; j++)
	sum2 += N[i][j];
      
      prob1 = (double) sum1 / (double) tot;
      prob2 = (double) sum2 / (double) tot;

      if (prob1 != 0)
	H += prob1 * log(prob1);
    
      if (prob2 != 0)
	H += prob2 * log(prob2);
    }
  
  return -H;
}


/* HPQ2
 * Returns the entropy induced by the points on the 
 * partition <c_s, c_t> and Q.
 * See line 13 of Algorithm 2 in SOM.
 *
 * Parameters
 *   N : q x p matrix containing the number of points
 *       for each cell of the grid formed by P and Q
 *       partitions.
 *   np : p vector containing the number of points
 *        for each cell of the grid formed by P
 *   q : number of rows of N (number of partitions in Q)
 *   p : number of cols of N and length of np (number of partitions in P)
 *   s : s in c_s
 *   t : t in c_t
 */
double HPQ2(int **N, int *np, int q, int p, int s, int t)
{
  int i, j;
  int sum, tot;
  double prob, H;
  
  
  if (s == t)
    return 0.0;

  tot = 0;
  for (i=s; i<t; i++)
    tot += np[i];

  H = 0.0;
  for (i=0; i<q; i++)
    {
      sum = 0;
      for (j=s; j<t; j++)
	sum += N[i][j];
      
      prob = (double) sum / (double) tot;

      if (prob != 0)
	H += prob * log(prob);
    }
  
  return -H;
}


/* EquipartitionYAxis
 * Returns the map Q: D -> {0, ...,y-1}. See Algorithm 3 in SOM.
 * 
 * Parameters
 *   Dy (IN): y-data sorted in increasing order
 *   n (IN): length of Dy
 *   y (IN): an integer greater than 1
 *   Qm (OUT): the map Q. Qm must be a preallocated vector 
 *       of size n.
 * Return
 *   q : the real value of y. It can be < y.
 */
int EquipartitionYAxis(double *Dy, int n, int y, int *Qm)
{
  int i, j, s, h, curr;
  double rowsize;
  
  for (i=0; i<n; i++)
    Qm[i] = -1;

  i = 0;
  curr = 0;
  rowsize =  (double) n / (double) y;
  while (i < n)
    {
      s = 1;
      for (j=i+1; j<n; j++)
	if (Dy[i] == Dy[j])
	  s++;
	else
	  break;
      
      h = 0;
      for (j=0; j<n; j++)
	if (Qm[j] == curr)
	  h++;
      
      if ((h != 0) && (fabs(h+s-rowsize) >= fabs(h-rowsize)))
	{
	  curr++;
	  rowsize = (double) (n-i) / (double) (y-curr);
	}

      for (j=0; j<s; j++)
	Qm[i+j] = curr;
      
      i += s;
    }

  return curr + 1;
}


/* GetClumpsPartition
 * Returns the map P: D -> {0, ...,k-1}.
 * 
 * Parameters
 *   Dx (IN) : x-data sorted in increasing order
 *   n (IN) : length of Dx
 *   Qm (IN) : the map Q computed by EquipartitionYAxis sorted in increasing 
 *        order by Dx-values.
 *   Pm (OUT) : the map P. Pm must be a preallocated vector of size n. 
 * Return
 *   k : number of clumps in Pm.
 */
int GetClumpsPartition(double *Dx, int n, int *Qm, int *Pm)
{
  int i, j, s, c, flag;
  int *Qm_tilde;
  
 
  Qm_tilde = (int *) malloc (n * sizeof(int));
  for(i=0; i<n; i++)
    Qm_tilde[i] = Qm[i];

  i = 0;
  c = -1; 
  while (i < n)
    {
      s = 1;
      flag = 0;
      for (j=i+1; j<n; j++)
	if (Dx[i] == Dx[j])
	  {
	    s++;
	    if (Qm_tilde[i] != Qm_tilde[j])
	      flag = 1;
	  }
	else
	  break;
      
      if ((s > 1) && (flag == 1))
	{
	  for (j=0; j<s; j++)
	    Qm_tilde[i+j] = c;
	  c--;
	}
      
      i += s;
    }
  
  i = 0;
  Pm[0] = i;
  for (j=1; j<n; j++)
    {
      if (Qm_tilde[j] != Qm_tilde[j-1])
	i++;
      Pm[j] = i;
    }
  
  free(Qm_tilde);
  return i+1;
}


/* GetSuperclumpsPartition
 * Returns the map P: D -> {0, ...,k-1}.
 * 
 * Parameters
 *   Dx (IN) : x-data sorted in increasing order
 *   n (IN) : length of Dx
 *   Qm (IN) : the map Q computed by EquipartitionYAxis sorted 
 *       in increasing order by Dx-values.
 *   k_hat (IN) : maximum number of clumps 
 *       Pm (IN): the map P. Pm must be a preallocated vector 
 *       of size n. 
 * Return
 *   k : number of clumps in Pm.
 */
int GetSuperclumpsPartition(double *Dx, int n, int *Qm, int k_hat, int *Pm)
{
  int i, k, p;
  double *Dp;
  
  /* compute clumps */
  k = GetClumpsPartition(Dx, n, Qm, Pm);

  if (k > k_hat) /* superclumps */
    {
      Dp = (double *) malloc (n * sizeof(double));
      for (i=0; i<n; i++)
	Dp[i] = (double) Pm[i];
      p = EquipartitionYAxis(Dp, n, k_hat, Pm);
      free(Dp);
      return p;
    }
  else
    return k;
 }


/* ApproxOptimizeXAxis
 * Returns the map P: D -> {0, ...,k-1}. See Algorithm 2 in SOM.
 * 
 * Parameters
 *   Dx (IN) : x-data sorted in increasing order by Dx-values
 *   Dy (IN) : y-data sorted in increasing order by Dx-values
 *   n (IN) : length of Dx and Dy
 *   Qm (IN) : the map Q computed by EquipartitionYAxis sorted 
 *   in increasing order by Dx-values.
 *   q (IN) : number of clumps in Qm
 *   Pm (IN) : the map P computed by GetSuperclumpsPartition 
 *   sorted in increasing order by Dx-values.
 *   p (IN) : number of clumps in Pm
 *   x (IN) : grid size on x-values 
 *   I (OUT) : the normalized mutual information vector. It 
 *   will contain I_{k,2}, ..., I_{k, x}. I must be a 
 *   preallocated array of dimension x-1.
 */
void ApproxOptimizeXAxis(double *Dx, double *Dy, int n, int *Qm, int q, 
			 int *Pm, int p, int x, double *I)
{
  int i, j, s, t, l;
  int **N; /* contains the number of samples for each cell Q x P*/
  int *np; /* contains the number of samples for each cell P */
  int *c; /* contains c_1, ..., c_k */
  double **IM; /* mutual information matrix, I_t,l */
  double f, fmax, r1, r2;
  double hq;
  double **hpq2;


  /* if p==1 return I=0 */
  if (p == 1)
    {
      for (i=0; i<x-1; i++)
	I[i] = 0.0;
      return;
    }
  
  /* alloc the N matrix */
  N = (int **) malloc (q * sizeof(int *));
  for (i=0; i<q; i++)
    {
      N[i] = (int *) malloc (p * sizeof(int));
      for (j=0; j<p; j++)
	N[i][j] = 0;
    }
  
  /* alloc the np vector */
  np = (int *) malloc (p * sizeof(int));
  for (j=0; j<p; j++)
    np[j] = 0;

  /* fill N and np */
  for (i=0; i<n; i++)
    {
      N[Qm[i]][Pm[i]] += 1;
      np[Pm[i]] += 1;
    }

  /* compute c_1, ..., c_k */
  c = (int *) malloc (p * sizeof(int));
  c[0] = np[0];
  for (i=1; i<p; i++)
    c[i] = np[i] + c[i-1];

  /* alloc the IM matrix */
  IM = (double **) malloc ((p+1) * sizeof(double *));
  for (i=0; i<=p; i++)
    {
      IM[i] = (double *) malloc ((x+1) * sizeof(double));
      for (j=0; j<=x; j++)
	IM[i][j] = 0.0;
    }
  
  /* compute H(Q)*/
  hq = HQ(N, q, p, n);
  
  /* Find the optimal partitions of size 2 */
  /* Algorithm 2 in SOM, lines 4-8 */
  for (t=2; t<=p; t++)
    {
      fmax = -DBL_MAX;
      for (s=1; s<=t; s++)
	{
	  f = HP3(np, p, s, t) - HPQ3(N, np, q, p, s, t);
	  if (f > fmax)
	    {
	      IM[t][2] = hq + f;
	      fmax = f;
	    }
	}
    }
    
  /* precomputed H(<c_s, c_t>, Q) matrix */
  hpq2 = (double **) malloc ((p+1) * sizeof(double *));
  for (i=0; i<=p; i++)
    {
      hpq2[i] = (double *) malloc ((p+1) * sizeof(double));
      for (j=0; j<=p; j++)
	hpq2[i][j] = 0.0;
    }
  for (t=3; t<=p; t++)
    for (s=2; s<=t; s++)
      hpq2[s][t] = HPQ2(N, np, q, p, s, t);
    
  /* inductively build the rest of the table of optimal partitions */
  /* Algorithm 2 in SOM, lines 11-17 */
  for (l=3; l<=x; l++)
    for (t=l; t<=p; t++)
      {
	fmax = -DBL_MAX;
	for (s=l-1; s<=t; s++)
	  {
	    r1 = (double) c[s-1] / (double) c[t-1];
	    r2 = (double) (c[t-1] - c[s-1]) / (double) c[t-1];
	    f = (r1 * (IM[s][l-1] - hq)) - (r2 * hpq2[s][t]);
	    
	    if (f > fmax)
	      {
		IM[t][l] = hq + f;
		fmax = f;
	      }
	  }
      }

  /* Algorithm 2 in SOM, line 19 */
  if (x > p)
    {
      for (i=p+1; i<=x; i++)
	IM[p][i] = IM[p][p];
    }
  
  /* fill I */
  for (i=2; i<=x; i++)
    I[i-2] = IM[p][i] / MIN(log(i), log(q));
    
  /* free */
  for (i=0; i<q; i++)
    free(N[i]);
  free(N);
  free(np);
  free(c);
  for (i=0; i<=p; i++)
    {
      free(IM[i]);
      free(hpq2[i]);
    }
  free(IM);
  free(hpq2);

  return;
}

/****** START QUICKSORT ******/

void swap(double *x, int *idx, int i, int j)
{
  double x_t;
  int idx_t;
  
  x_t = x[i];
  x[i] = x[j];
  x[j] = x_t;

  idx_t = idx[i];
  idx[i] = idx[j];
  idx[j] = idx_t;
}

/* sort x from index l to index u and idx according to x */
void quicksort(double *x, int *idx, int l, int u)
{
  int i, m;

  if (l >= u)
    return;

  m = l;
  for (i=l+1; i<=u; i++)
    if (x[i] < x[l])
      swap(x, idx, ++m, i);
  swap(x, idx, l, m);
  quicksort(x, idx, l, m-1);
  quicksort(x, idx, m+1, u);
}

/* sort
 * Sort x and idx (of length n) according to x.
 */
void sort(double *x, int *idx, int n)
{
  quicksort(x, idx, 0, n-1);
}

/****** END QUICKSORT ******/
