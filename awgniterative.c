/*
  	FILENAME: awgniterative.c
  	AUTHOR: Tadashi Wadayama

	NAME: awgniterative

	SYNOPSYS: awgniterative spmat_file snr (maxi seed stop #err disp)

	Output format of simulation results
	snr pb pB var #ebits #bits #eblks #blks aveitr 
	seed maxi N M file stop #err

	Symbols

	snr   : Eb/N0
	pb    : bit error probability after decoding 
	pB    : block error probability after decoding 
	var   : variance of the noise
	#ebits: number of error bits
	#bits : number of total bits generated
	#eblks: number of error blocks
	#blks : number of total blocks generated
	aveitr: average number of iterations per block
	seed  : the seed of the random number generator
	maxi  : the maximum number of iteration 
	N     : code length
	M     : number of parity bits
	file  : file name of parity check matrix (in spmat form)
	stop  : 
	        0->simulation stops when #ebits becomes #err
	        1->simulation stops when #eblks becomes #err
	#err  : number of errors enough to stop a simulation
	disp  : display mode(disp = 1: display)
	
	Caution: The error bits are counted over a whole word.
	It may not be a traditional definition of bit error 
	probability after decoding.

	DESCRIPTION:
	
	The program is a simulation program 
	for iterative decoding for low density 
	parity check codes. An additive white Gaussian noise channel 
	is assumed.

	EXAMPLE:
	awgniterative 981.500 2.0 (10 123 0 100 0)

	HOW TO MAKE:
	gcc -O2 -o awgniterative awgniterative.c -lm

	HISTORY:
  	SINCE : Jan. 17, 2000
	Jan.19: bug fix about malloc of bcjr_for0,etc.

	Copyright (C) Tadashi Wadayama

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* =============================================== */
/*  Each item corresponds to bit with value 1 in H */
/* =============================================== */

typedef struct _ITEM
{
  int n,m;			/* (m,n):the position in H */
  void *right;		/* The pointer to next item in the same row */
  void *down;		/* The pointer to next item in the same column*/
  double r0,r1;
  double q0,q1;
} item;

/* =============================================== */
/*  The following is the central data structure.   */
/*  It represnts a sparce matrix by using linked   */
/*  list.  Every item in the lists have            */
/*  the following links to each other              */
/*

 item -> item -> ...
  |       |
  V       V
 item -> item -> ...
  |       |
  V       V
/* =============================================== */

typedef struct _SPMATRIX
{
  int N;			/* number of colmun */
  int M;			/* number of row */
  int* num_ones_in_col;		/* number of ones in a column of H */
  int* num_ones_in_row;		/* number of ones in a row of H */
  int biggest_num_ones_col;	/* The biggest number of ones in column */
  int biggest_num_ones_row;	/* The biggest number of ones in row */
  item* start_col_list;		/* a column list starts from here*/
  item* start_row_list;		/* a row list starts from here*/
} sparce_matrix;

/* =============================================== */
/*  The following data structure is used for       */
/*  storing information related to the simulation. */
/* =============================================== */

typedef struct _SIM_PARA
{
  int N;			/* code length */
  int M;			/* the number of redundancy */
  double snr;			/* Eb/N0 */
  double var;			/* variance of the noise */
  double* rword;		/* received_word */

				/* for BCJR algorithm */
  double* bcjr_for0;		/* forward probabilty of 0-state */
  double* bcjr_for1;		/* forward probabilty of 1-state */
  double* bcjr_back0;		/* backword probabilty of 0-state */
  double* bcjr_back1;		/* backword probabilty of 1-state */
  double* bcjr_tmp0;		/* likelihood for symbol 0 */
  double* bcjr_tmp1;		/* likelihood for symbol 1 */
  double* bcjr_tmp_q0;
  double* bcjr_tmp_q1;
				/* for updown algorithm */
  double* ud_upward0;		/* upward probabilty of 0-state */
  double* ud_upward1;		/* upward probabilty of 1-state */
  double* ud_downward0;		/* downward probabilty of 0-state */
  double* ud_downward1;		/* downward probabilty of 1-state */
  double* ud_tmp0;
  double* ud_tmp1;

  double* tmp_q0;		/* pseudo probability for 0 */
  double* tmp_q1;		/* pseudo probability for 1 */
  int* tmp_decision;		/* temporary decision for each bit */

  int max_iteration;		/* maximum number of iterations */
  int total_blocks;		/* number of transmitted blocks */
  int error_blocks;		/* number of error blocks */
  int total_bits;		/* number of transmitted bits */
  int error_bits;		/* number of error bits */
  int error_weight_in_word;	/* number of errors within a block */
  int seed;			/* seed for the random number generator */
  int num_iteration;		/* total number of required iterations */
  int stop;			/* stop flag*/
  int stop_err;			/* number of errors enough to stop */
  int display;			/* 0->non display mode, 1->display mode */
} simulation_parameters;

/* ================================================== */
/*  A normal Gaussian noise generator                 */
/* ================================================== */

double nrnd(double var)
{
  static int sw = 0;
  static double r1, r2, s;

  if (sw == 0) {
    sw = 1;
    do {

      r1 = 2 * drand48() - 1;
      r2 = 2 * drand48() - 1;

      s = r1 * r1 + r2 * r2;
    } while (s > 1 || s == 0);
    s = sqrt(-2 * log(s) / s);
    return r1 * s * sqrt(var);
  } else {
    sw = 0;
    return r2 * s * sqrt(var);
  }
}


/* ================================================== */
/*  The function shows the contents in a sparce_matrix*/
/* ================================================== */

void print_parameter(sparce_matrix* s)
{
  int m,i;
  item* p;
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      printf("\t(%d,%d:r0=%f,r1=%f) ",
	     p->m, p->n,p->r0,p->r1);
      p = p->right;
    }
    printf("\n");

    p = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      printf("\t(%d,%d:q0=%f,q1=%f) ",
	     p->m, p->n,p->q0,p->q1);
      p = p->right;
    }
    printf("\n\n");
  }
}

/* ================================================== */
/*  The function reads spmat_file and sets the        */
/*  sparce_matrix data structure.                     */
/*
This is a exmple of an spmat format.

16 12	                          N,M
4 3                               biggest_num_ones_row biggest_num_ones_col
4 4 4 4 4 4 4 4 4 4 4 4           num_ones_in_row
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3   num_ones_in_col
3 8 10 13                         row form (caution:colum no started from 1!)
4 7 9 13
2 5 7 10
4 6 11 14
3 9 15 16
1 6 9 10
4 8 12 15
2 6 12 16
1 7 14 16
3 5 12 14
2 11 13 15
1 5 8 11                    

The above spmat format corresponds to the following matrix:

parity_check_matrix
16
4
 0 0 1 0 0 0 0 1 0 1 0 0 1 0 0 0
 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0
 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0
 0 0 0 1 0 1 0 0 0 0 1 0 0 1 0 0
 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1
 1 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0
 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0
 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1
 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1
 0 0 1 0 1 0 0 0 0 0 0 1 0 1 0 0
 0 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0
 1 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0
*/
/* ================================================== */

void read_spmat_file(FILE* fp, sparce_matrix* a)
{
  int i,j;
  int n,m;
  char buf[128];
  item *new_item, *last_item;
  item* p;
  item* q;
  int tmp;
    
  fscanf(fp,"%d %d\n",&(a->N),&(a->M));	/* reading N and M */
  fscanf(fp,"%d %d\n",&(a->biggest_num_ones_row),&(a->biggest_num_ones_col));

				/* Initialization  */
  
  if ((a->num_ones_in_row = (int*)malloc(sizeof(int)*(a->M))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((a->num_ones_in_col = (int*)malloc(sizeof(int)*(a->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((a->start_col_list = (item*)malloc(sizeof(item)*(a->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((a->start_row_list = (item*)malloc(sizeof(item)*(a->M))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  for (i = 0; i <= a->M-1; i++) {
    fscanf(fp,"%d",&(a->num_ones_in_row[i]));
  }
  for (i = 0; i <= a->N-1; i++) {
    fscanf(fp,"%d",&(a->num_ones_in_col[i]));
    a->start_col_list[i].down = NULL; /* This line is appended.  */
  }

				/* making row links */
  for (i = 0; i <= a->M-1; i++) {
				/* processing i-th row */
    last_item = &(a->start_row_list[i]);
    for (j = 0; j <= a->num_ones_in_row[i]-1; j++) {
      fscanf(fp,"%d",&tmp);
      tmp--;
				/* making new item */
      if ((new_item = (item*)malloc(sizeof(item))) == NULL) {
	fprintf(stderr,"Can't allocate memory\n");
	exit(-1);
      }
      new_item->m = i;
      new_item->n = tmp;
      new_item->right = NULL;
      new_item->down = NULL;
      last_item->right = new_item;
      last_item = new_item;
    }
  }
				/* making column links */
  for (m = 0; m <= a->M-1; m++) {
    p = (a->start_row_list[m]).right;
    for (i = 0; i <= a->num_ones_in_row[m]-1; i++) {
      q = &(a->start_col_list[p->n]);
      while(q->down != NULL) {
	q = q->down;
      }
      q->down = p; 
      p = p->right;
    }
  }
}

/* ================================================== */
/*  AWGN channel                                      */
/* ================================================== */

void awgn_channel(sparce_matrix* s, simulation_parameters* param)
{
  int i;
  for (i = 0; i <= param->N-1; i++) {
    param->rword[i] = 1.0 + nrnd(param->var);
  }
}

/* ================================================== */
/*  The BCJR algorithm                                */
/* ================================================== */

void bcjr(sparce_matrix* s, simulation_parameters* p,int row_num)
{
  int i,m;
  item* pos;
  double t0,t1,sum;

				/* copy */
  /*

    likelihood function 
    exp(y x/ var)  for x = +1,-1
    y: received symbol
    x: code symbol

  */

  pos = (s->start_row_list[row_num]).right;
  for (i = 0; i <= s->num_ones_in_row[row_num]-1; i++) {
    p->bcjr_tmp0[i] = exp(p->rword[pos->n] / p->var);
    p->bcjr_tmp1[i] = exp(-p->rword[pos->n] / p->var);
    p->bcjr_tmp_q0[i] = pos->q0;
    p->bcjr_tmp_q1[i] = pos->q1;
    pos = pos->right;
  }
				/* forward computation */
  p->bcjr_for0[0] = 1.0;
  p->bcjr_for1[0] = 0.0;

  for (i = 1; i <= s->num_ones_in_row[row_num]; i++) {
    t0  =  p->bcjr_tmp0[i-1] * p->bcjr_tmp_q0[i-1] * p->bcjr_for0[i-1]
      + p->bcjr_tmp1[i-1] * p->bcjr_tmp_q1[i-1] * p->bcjr_for1[i-1];

    t1  = p->bcjr_tmp1[i-1] * p->bcjr_tmp_q1[i-1] * p->bcjr_for0[i-1]
      + p->bcjr_tmp0[i-1] * p->bcjr_tmp_q0[i-1] * p->bcjr_for1[i-1];

				/* scaling */
    sum = t0 + t1;
    p->bcjr_for0[i] = t0/sum;
    p->bcjr_for1[i] = t1/sum;
  }
  p->bcjr_for1[s->num_ones_in_row[row_num]] = 0.0;

				/* backward computation */

  p->bcjr_back0[s->num_ones_in_row[row_num]] = 1.0;
  p->bcjr_back1[s->num_ones_in_row[row_num]] = 0.0;

  for (i = s->num_ones_in_row[row_num]-1; i >= 0; i--) {
    t0 = p->bcjr_tmp0[i] * p->bcjr_tmp_q0[i] * p->bcjr_back0[i+1]
      + p->bcjr_tmp1[i] * p->bcjr_tmp_q1[i] * p->bcjr_back1[i+1];

    t1 = p->bcjr_tmp1[i] * p->bcjr_tmp_q1[i] * p->bcjr_back0[i+1]
      + p->bcjr_tmp0[i] * p->bcjr_tmp_q0[i] * p->bcjr_back1[i+1];

				/* scaling */
    sum = t0 + t1;
    p->bcjr_back0[i] = t0/sum;
    p->bcjr_back1[i] = t1/sum;
  }
  p->bcjr_back1[0] = 0.0;

				/* update for r0 and r1 */

  pos = (s->start_row_list[row_num]).right;
  for (i = 0; i <= s->num_ones_in_row[row_num]-1; i++) {

				/* extrinsic values */
    pos->r0 = p->bcjr_for0[i] * p->bcjr_back0[i+1]
      + p->bcjr_for1[i] * p->bcjr_back1[i+1];
    
    pos->r1= p->bcjr_for0[i] * p->bcjr_back1[i+1]
      + p->bcjr_for1[i] * p->bcjr_back0[i+1];
    pos = pos->right;
  }
}

/* ================================================== */
/*  Updown algorith for updating q0 and q1            */
/* ================================================== */

void updown(sparce_matrix* s, simulation_parameters* p,int col_num)
{
  int i,n;
  item* pos;
  double t0, t1;
  double sum;
				/* copy */
  pos = (s->start_col_list[col_num]).down;
  for (i = 0; i <= s->num_ones_in_col[col_num]-1; i++) {
    p->ud_tmp0[i] = pos->r0;
    p->ud_tmp1[i] = pos->r1;
    pos = pos->down;
  }
				/* downward computation */
  p->ud_downward0[0] = 1.0;
  p->ud_downward1[0] = 1.0;
  for (i = 1; i <= s->num_ones_in_col[col_num]; i++) {
    t0 = p->ud_tmp0[i-1] * p->ud_downward0[i-1];
    t1 = p->ud_tmp1[i-1] * p->ud_downward1[i-1];
    sum = t0 + t1;    
    p->ud_downward0[i] = t0/sum;
    p->ud_downward1[i] = t1/sum;
  }
				/* upward computation */
  p->ud_upward0[s->num_ones_in_col[col_num]] = 1.0;
  p->ud_upward1[s->num_ones_in_col[col_num]] = 1.0;
  for (i = s->num_ones_in_col[col_num]-1; i >= 0; i--) {
    t0 = p->ud_tmp0[i] * p->ud_upward0[i+1];
    t1 = p->ud_tmp1[i] * p->ud_upward1[i+1];
    sum = t0 + t1;
    p->ud_upward0[i] = t0/sum;
    p->ud_upward1[i] = t1/sum;
  }

				/* update for q0 and q1 */
  pos = (s->start_col_list[col_num]).down;
  for (i = 0; i <= s->num_ones_in_col[col_num]-1; i++) {
    t0 = p->ud_downward0[i] * p->ud_upward0[i+1];
    t1 = p->ud_downward1[i] * p->ud_upward1[i+1];

				/* scaling */
    sum = t0 + t1;
    pos->q0 = t0/sum;
    pos->q1 = t1/sum;

    if (pos->q0 == 0) {
      //      printf("caution t0  %16.12f\n",pos->q0);
      pos->q0 = 1E-8;
      pos->q1 = 1.0 - 1E-8;
    }
    if (pos->q1 == 0) {
      //      printf("caution t1  %16.12f\n",pos->q1);
      pos->q1 = 1E-8;
      pos->q0 = 1.0 - 1E-8;
    }

    pos = pos->down;
  }
				/* update pseudo probability for each bit */
  
  t0 = exp(p->rword[col_num] / p->var) * p->ud_upward0[0];
  t1 = exp(-p->rword[col_num] / p->var) * p->ud_upward1[0];

				/* scaling */

  sum = t0 + t1;
  p->tmp_q0[col_num] = t0/sum;
  p->tmp_q1[col_num] = t1/sum;

				/* temporary decision for each bit*/


  if (p->tmp_q1[col_num] > p->tmp_q0[col_num]) 
    p->tmp_decision[col_num] = 1;
  else 
    p->tmp_decision[col_num] = 0;

				/* error count */
  p->error_weight_in_word += p->tmp_decision[col_num];
}

/* ================================================== */
/*  parity check function                             */
/*  return value = 0 : tmp_decision is a codeword     */
/*  return value = 1 : tmp_decision is not a codeword */
/* ================================================== */

int parity_check(sparce_matrix* s, simulation_parameters* p)
{
  int m,i;
  item* pos;
  int parity;

  for (m = 0; m <= s->M-1; m++) {
    parity = 0;
    pos = (s->start_row_list[m]).right;
    for (i = 0; i <= s->num_ones_in_row[m]-1; i++) {
      parity = (parity + p->tmp_decision[pos->n]) % 2;
      pos = pos->right;
    }
    if (parity == 1) return 1;
  }
  return 0;
}

/* ================================================== */
/*  print function for simulation results             */
/* ================================================== */

void print_results(FILE* out,char* file, simulation_parameters* param)
{
				/* print simulation results */
  if (param->display == 0)
    fprintf(out,
	    "#snr pb pB var #eblks #blks #ebits #bits aveitr seed maxitr\n");
  fprintf(out,"%16.12e %16.12e %16.12e %16.12e %d %d %d %d %f %d %d %d %d %s %d %d\n",
	  param->snr,
	  (double)param->error_bits/param->total_bits,
	  (double)param->error_blocks/param->total_blocks,
	  param->var,
	  param->error_bits,
	  param->total_bits,
	  param->error_blocks,
	  param->total_blocks,
	  (double)param->num_iteration/param->total_blocks,
	  param->seed,
	  param->max_iteration,
	  param->N,
	  param->M,
	  file,
	  param->stop,
	  param->stop_err
	  );
}

/* ================================================== */
/*  sum product decoder                               */
/*  return value = 0: success                         */
/*  return value = 1: failure                         */
/* ================================================== */

int sum_product_decoder(sparce_matrix* s, simulation_parameters* param)
{
  int i,n,m;
  item* pos;
  int parity;
				/* initialize q0 and q1 to be 1.0 */
  for (m = 0; m <= s->M-1; m++) {
    pos = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      pos->q0 = 1.0;
      pos->q1 = 1.0;
      pos = pos->right;
    }
  }
				/* iterative decoding */
  for (i = 1; i <= param->max_iteration; i++) {
    param->error_weight_in_word = 0;
    param->num_iteration++;
				/* row processing */

    for (m = 0; m <= s->M-1; m++) bcjr(s,param,m); 

				/* column processing */
    for (n = 0; n <= s->N-1; n++) updown(s,param,n);

				/* parity check for temporary decision */
    parity = parity_check(s,param); 

    if (parity == 0) return 0;	/* successful decoding */
  }
  return 1;			/* decoding failure */
}

/* ================================================== */
/*  Initialization on simulation parameters           */
/* ================================================== */

void init_simulation_param(sparce_matrix* s, simulation_parameters* param)
{
  param->N = s->N;
  param->M = s->M;
  param->total_blocks = 0;
  param->error_blocks = 0;
  param->total_bits = 0;
  param->error_bits = 0;
  param->num_iteration = 0;

  if ((param->rword = (double*)malloc(sizeof(double)*param->N)) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->bcjr_for0 
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_row + 1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
  if ((param->bcjr_for1 
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_row + 1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
  if ((param->bcjr_back0
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_row + 1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
  if ((param->bcjr_back1
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_row + 1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->bcjr_tmp0
       = (double*)malloc(sizeof(double)*s->biggest_num_ones_row)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->bcjr_tmp1
       = (double*)malloc(sizeof(double)*s->biggest_num_ones_row)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->bcjr_tmp_q0
       = (double*)malloc(sizeof(double)*s->biggest_num_ones_row)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->bcjr_tmp_q1
       = (double*)malloc(sizeof(double)*s->biggest_num_ones_row)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->ud_upward0
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_col+1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->ud_upward1
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_col+1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->ud_downward0
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_col+1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->ud_downward1
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_col+1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->ud_tmp0
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_col+1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->ud_tmp1
       = (double*)malloc(sizeof(double)*(s->biggest_num_ones_col+1))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->tmp_q0
       = (double*)malloc(sizeof(double)*s->N)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->tmp_q1
       = (double*)malloc(sizeof(double)*s->N)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((param->tmp_decision
       = (int*)malloc(sizeof(int)*param->N)) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
}

/* ================================================== */
/*  Count error blocks and bits                       */
/* ================================================== */

void error_count(simulation_parameters* param)
{
  int i;

  if (param->error_weight_in_word != 0) 
    param->error_blocks++;
  param->error_bits += param->error_weight_in_word;
}


/* ================================================== */
/*  The function shows contents of a sparce_matrix    */
/*  data in spmat form.                               */
/* ================================================== */

void print_spmatrix_in_spmatform(sparce_matrix* s)
{

  /*

This is a exmple of an spmat format.

16 12	                          N,M
4 3                               biggest_num_ones_row biggest_num_ones_col
4 4 4 4 4 4 4 4 4 4 4 4           num_ones_in_row
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3   num_ones_in_col
3 8 10 13                         row form
4 7 9 13
2 5 7 10
4 6 11 14
3 9 15 16
1 6 9 10
4 8 12 15
2 6 12 16
1 7 14 16
3 5 12 14
2 11 13 15
1 5 8 11                    

The above spmat format corresponds to the following matrix:

parity_check_matrix
16
4
 0 0 1 0 0 0 0 1 0 1 0 0 1 0 0 0
 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0
 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0
 0 0 0 1 0 1 0 0 0 0 1 0 0 1 0 0
 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1
 1 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0
 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0
 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1
 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1
 0 0 1 0 1 0 0 0 0 0 0 1 0 1 0 0
 0 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0
 1 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0

   */

  int m,i;
  item* p;
  printf("%d %d\n",s->N,s->M);
  printf("%d %d\n",s->biggest_num_ones_row,s->biggest_num_ones_col);
  for (i = 0; i <= s->M-1; i++) {
    printf("%d ",s->num_ones_in_row[i]);
  }
  printf("\n");
  for (i = 0; i <= s->N-1; i++) {
    printf("%d ",s->num_ones_in_col[i]);
  }
  printf("\n");
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      printf("%d ",(p->n)+1);
      p = p->right;
    }
    printf("\n");
  }
  
}

void print_spmatrix_in_colform(sparce_matrix* s)
{
  int n,i;
  item* p;

  printf("N = %d\n",s->N);
  printf("M = %d\n",s->M);

  printf("num_ones_in_row\n");
  for (i = 0; i <= s->M-1; i++) {
    printf("%d ",s->num_ones_in_row[i]);
  }
  printf("\n");

  printf("num_ones_in_col\n");
  for (i = 0; i <= s->N-1; i++) {
    printf("%d ",s->num_ones_in_col[i]);
  }
  printf("\n");
  printf("biggest_num_ones_row = %d\n",s->biggest_num_ones_row);
  printf("biggest_num_ones_col = %d\n",s->biggest_num_ones_col);

  printf("COLUMN form\n");
  for (n = 0; n <= s->N-1; n++) {
    p = (s->start_col_list[n]).down;
    printf("%d-th col: ",n);
    for (i = 1; i <= s->num_ones_in_col[n]; i++) {
      printf("(%d,%d)->",p->m, p->n);
      p = p->down;
    }
    printf("\n");
  }
}

int main(int argc,char **argv)
{
  FILE* fp;
  sparce_matrix s;
  simulation_parameters param;
  int i;
				/* user interface */
  if (argc < 3) {
    printf("usage : awgniterative file snr (maxitr seed stop #err disp)\n");
    printf("file  : parity check matrix (in spmat format)\n");
    printf("snr   : Eb/N0\n");
    printf("maxitr: maximum number of iterations\n");
    printf("seed  : seed for the random number generator\n");
    printf("stop  : = 0: simulation stops when #ebits becomes #err\n");
    printf("        = 1: simulation stops when #eblks becomes #err\n");
    printf("disp  : = 0: non display\n");
    printf("        = 1: display\n");
    exit(-1);
  }
  if ((fp = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",argv[1]);
    exit(-1);
  }
  param.snr = atof(argv[2]);
  if (argc >= 4) param.max_iteration = atoi(argv[3]);
  else param.max_iteration = 20; /* default maxi */

  if (argc >= 5) param.seed = atoi(argv[4]);
  else param.seed = 1234;	/* default seed */

  if (argc >= 6) param.stop = atoi(argv[5]);
  else param.stop = 0;		/* default:bit error criteria */

  if (argc >= 7) param.stop_err = atoi(argv[6]);
  else param.stop_err = 100;	/* 100 errors  */

  if (argc >= 8) param.display = atoi(argv[7]);
  else param.display = 0;	/* default: non-display */

				/* initialize random number generator */
  srand48(param.seed);
				/* reading spmat file */
  read_spmat_file(fp, &s);
  /* print_spmatrix_in_spmatform(&s); */

  /*
    Assumption:
    * binary symbol 0 is mapped to +1.
    * binary symbol 1 is mapped to -1.
    * transmitted word is (+1,+1,....,+1).
    * The variance of noise is param.var.
   */
				/* compute variance */

  param.var = 0.5 * (1.0/pow(10.0,param.snr/10.0)) 
    * (double)s.N/(double)(s.N-s.M);
				/* initialization of simulation parameters */
  init_simulation_param(&s,&param);

				/* simulation loop */
  while(1) {
    param.total_blocks++;
    param.total_bits += param.N;
    
    awgn_channel(&s,&param);
    sum_product_decoder(&s,&param);
    error_count(&param);
    if (param.display == 1) {
      print_results(stderr,"****",&param);
    }
				/* stop criteria 0 (stop = 0) */
    if ((param.stop == 0) && (param.error_bits >= param.stop_err)) break;

				/* stop criteria 1 (stop = 1) */
    if ((param.stop == 1) && (param.error_blocks >= param.stop_err)) break;
  }
				/* print results */
  print_results(stdout,argv[1],&param);

}




