/*
  	FILENAME: mkencoder.c
  	AUTHOR: Tadashi Wadayama

	NAME:
	mkencoder

	SYNOPSYS: mkencoder spmat_file 
	
	spmat_file : a parity check matrix in spmat format

	DESCRIPTION:
	
	The program makes encoder matrix from a 
	a given parity_check_matrix in spmat form.

	EXAMPLE:
	mkencoder mkencoder.test

	BUGS:
	

	RELATED DOCUMENTS or SOURCE CODES:
	../spmat.form

	HOW TO MAKE:
	gcc -O2 -o mkencoder mkencoder.c

	HISTORY:
  	SINCE : Jan. 31, 2000
	Feb.8  bug fix : num_ones_in_col in decoder matrix
	Feb.11 bug fix : biggest_num_ones_row in encoder matrix
*/

#include <stdio.h>
#include <string.h>

int DISPLAY;

/* =============================================== */
/*  Each item corresponds to bit with value 1 in H */
/* =============================================== */

typedef struct _ITEM
{
  int n,m;			/* (m,n):the position in H */
  void *right;		/* The pointer to next item in the same row */
  void *down;		/* The pointer to next item in the same column*/
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

*/
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

/* ================================================== */
/*  The function shows contents of a sparce_matrix    */
/*  data in row form.                                 */
/*  The function can be used for themplete for row    */
/*  traverse.                                         */
/* ================================================== */

void print_spmatrix_in_rowform(sparce_matrix* s)
{
  int m,i;
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

  printf("ROW form\n");
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    printf("%d-th row: ",m);
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      printf("(%d,%d)->",p->m, p->n);
      p = p->right;
    }
    printf("\n");
  }
}

/* ================================================== */
/*  The function shows contents of a sparce_matrix    */
/*  data in column form.                              */
/*  The function can be used for themplete for column */
/*  traverse.                                         */
/* ================================================== */

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

/* ================================================== */
/*  The function shows contents of a sparce_matrix    */
/*  data in spmat form.                               */
/* ================================================== */

void print_spmatrix_in_spmatform(FILE* fp, sparce_matrix* s)
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
  fprintf(fp,"%d %d\n",s->N,s->M);
  fprintf(fp,"%d %d\n",s->biggest_num_ones_row,s->biggest_num_ones_col);
  for (i = 0; i <= s->M-1; i++) {
    fprintf(fp,"%d ",s->num_ones_in_row[i]);
  }
  fprintf(fp,"\n");
  for (i = 0; i <= s->N-1; i++) {
    fprintf(fp,"%d ",s->num_ones_in_col[i]);
  }
  fprintf(fp,"\n");
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      fprintf(fp,"%d ",(p->n)+1);
      p = p->right;
    }
    fprintf(fp,"\n");
  }
  
}


void print_result(sparce_matrix* s, sparce_matrix* t, 
		  char* encoder, char* decoder)
{
  int n,m,i;
  item* p;
  int* col_perm;
  int* inv_perm;
  FILE* enc;
  FILE* dec;
  int tmp;

  if ((enc = fopen(encoder,"w")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",encoder);
    exit(-1);
  }

  if ((dec = fopen(decoder,"w")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",decoder);
    exit(-1);
  }

  if ((col_perm = (int*)malloc(sizeof(int)*(s->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((inv_perm = (int*)malloc(sizeof(int)*(s->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  for (i = 0; i <= s->N-1; i++) col_perm[i] = i;

  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;

    if (p->n != m) {
      tmp = col_perm[m];
      col_perm[m] = col_perm[p->n];
      col_perm[p->n] = tmp;
    }
  }
  for (n = 0; n <= s->N-1; n++) {
    inv_perm[col_perm[n]] = n;
  }


  for (n = 0; n <= s->N-1; n++) s->num_ones_in_col[n] = 0;
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      s->num_ones_in_col[inv_perm[p->n]]++;
      p = p->right; 
    }
  }

  s->biggest_num_ones_col = 0;
  for (n = 0; n <= s->N-1; n++) {
    if (s->biggest_num_ones_col < s->num_ones_in_col[n])
      s->biggest_num_ones_col = s->num_ones_in_col[n];
  }

  s->biggest_num_ones_row = 0;
  for (m = 0; m <= s->M-1; m++) {
    if (s->biggest_num_ones_row < s->num_ones_in_row[m])
      s->biggest_num_ones_row = s->num_ones_in_row[m];
  }
				/* printing encoder matrix  */

  fprintf(enc,"%d %d\n",s->N,s->M);
  fprintf(enc,"%d %d\n",s->biggest_num_ones_row,s->biggest_num_ones_col);
  for (i = 0; i <= s->M-1; i++) {
    fprintf(enc,"%d ",s->num_ones_in_row[i]);
  }
  fprintf(enc,"\n");
  for (i = 0; i <= s->N-1; i++) {
    fprintf(enc,"%d ",s->num_ones_in_col[i]);
  }
  fprintf(enc,"\n");
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    for (i = 1; i <= s->num_ones_in_row[m]; i++) {
      fprintf(enc,"%d ",inv_perm[p->n]+1);  
      p = p->right;
    }
    fprintf(enc,"\n");
  }

  for (n = 0; n <= t->N-1; n++) t->num_ones_in_col[n] = 0;
  for (m = 0; m <= t->M-1; m++) {
    p = (t->start_row_list[m]).right;
    for (i = 1; i <= t->num_ones_in_row[m]; i++) {
      t->num_ones_in_col[inv_perm[p->n]]++;
      p = p->right;
    }
  }
				  /* printing decoder matrix  */
  fprintf(dec,"%d %d\n",t->N,t->M);
  fprintf(dec,"%d %d\n",t->biggest_num_ones_row,t->biggest_num_ones_col);
  for (i = 0; i <= t->M-1; i++) {
    fprintf(dec,"%d ",t->num_ones_in_row[i]);
  }
  fprintf(dec,"\n");
  for (i = 0; i <= t->N-1; i++) {
    fprintf(dec,"%d ",t->num_ones_in_col[i]);  
  }
  fprintf(dec,"\n");
  for (m = 0; m <= t->M-1; m++) {
    p = (t->start_row_list[m]).right;
    for (i = 1; i <= t->num_ones_in_row[m]; i++) {
      fprintf(dec,"%d ",inv_perm[p->n]+1);  
      p = p->right;
    }
    fprintf(dec,"\n");
  }

}

/* ================================================== */
/*  The function reads spmat_file and sets the        */
/*  sparce_matrix data structure.                     */
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
    for (i = 1; i <= a->num_ones_in_row[m]; i++) {
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
/*  Gaussian ellimination                             */
/* ================================================== */

void gaussian(FILE* fp, sparce_matrix* a)
{
  int i,j,k;
  int n,m;
  char buf[128];
  item *new_item, *last_item;
  item* p;
  item* q;
  int tmp;
  int* current_row;
  int w;
  int leader;
  int eleader;

 
  fscanf(fp,"%d %d\n",&(a->N),&(a->M));	/* reading N and M */
  fscanf(fp,"%d %d\n",&(a->biggest_num_ones_row),&(a->biggest_num_ones_col));

				/* Initialization  */

  if ((current_row = (int*)malloc(sizeof(int)*(a->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
  
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
  }

  for (i = 0; i <= a->M-1; i++) {

                                         /* i-th row process */
    if (DISPLAY == 1) printf("now processing %d-th row\n",i);
    for (j = 0; j <= a->N-1; j++) current_row[j] = 0;
    for (j = 0; j <= a->num_ones_in_row[i]-1; j++) {
      fscanf(fp,"%d",&tmp);
      current_row[tmp-1] = 1;
    }
  retry:
				/* print row */
    if (DISPLAY == 1) {
      for (j = 0; j <= a->N-1; j++) printf("%d ",current_row[j]);printf("\n");
    }

	                       /* find leader  */
    k = 0; 
    while ((current_row[k] == 0) && (k <= a->N-1)) k++; 
    leader = k;
    if (DISPLAY == 1) printf("current_row leader = %d\n",leader);

	                       /* ellimination loop */
    for(j = 0; j <= i-1; j++) {
      if (DISPLAY == 1) printf("\t%d-th row: ",j);

      p = (a->start_row_list[j]).right;
      if (p == NULL) eleader = -1;
      else eleader = p->n;

      if (DISPLAY == 1) printf("above  leader      = %d\n",eleader);

      if (leader == eleader) {
	if (DISPLAY == 1) printf("elimination process\n");
	for (k = 1; k <= a->num_ones_in_row[j]; k++) {
	  current_row[p->n] = (current_row[p->n]+1) % 2;
	  p = p->right;
	}
	goto retry;
      }
    }
    
                                       /* count row weight */
    w = 0;
    for (j = 0; j <= a->N-1; j++) w += current_row[j];
    a->num_ones_in_row[i] = w;
    if (w == 0) {
      fprintf(stderr,"The input matrix is not full rank!\n");
      exit(-1);
    }
    if (DISPLAY == 1) {
      printf("***** final result of %d-th row***************************\n",i);
      for (j = 0; j <= a->N-1; j++) printf("%d ",current_row[j]);
      printf("w = %d\n",w);
      printf("***********************************************************\n");
    }

				/* making row links */
    last_item = &(a->start_row_list[i]);
    k = 0;
    for (j = 0; j <= a->num_ones_in_row[i]-1; j++) {
      while ((current_row[k] == 0) && (k <= a->N-1)) k++;
				/* making new item */
      if ((new_item = (item*)malloc(sizeof(item))) == NULL) {
	fprintf(stderr,"Can't allocate memory\n");
	exit(-1);
      }
      new_item->m = i;
      new_item->n = k;
      new_item->right = NULL;
      new_item->down = NULL;
      last_item->right = new_item;
      last_item = new_item;
      k++;
    }
  }
}

int row_perm(sparce_matrix* s)
{
  int m,i,j;
  item *p;
  int* perm_order;
  int tmp;
  item itmp;

  if ((perm_order = (int*)malloc(sizeof(int)*s->M)) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
  for (m = 0; m <= s->M-1; m++) {
    p = (s->start_row_list[m]).right;
    perm_order[m] = p->n;
  }

  for (i = 0; i <= s->M-1; i++) {
    for (j = i+1; j <= s->M-1; j++) {
      if (perm_order[i] > perm_order[j]) {
				/* swap */
	tmp = perm_order[i];
	perm_order[i] = perm_order[j];
	perm_order[j] = tmp;
	
	tmp = s->num_ones_in_row[i];
	s->num_ones_in_row[i] = s->num_ones_in_row[j];
	s->num_ones_in_row[j] = tmp;

	itmp = s->start_row_list[i];
	s->start_row_list[i] = s->start_row_list[j];
	s->start_row_list[j] = itmp;
      }
    }
  }

}


int main(int argc,char **argv)
{
  FILE* fp;
  sparce_matrix s,t;
  char encoder[128];
  char decoder[128];
  int i;

  if (argc < 2) {
    printf("usage: mkencoder spmat_file\n");
    exit(-1);
  }
  if ((fp = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",argv[1]);
    exit(-1);
  }

  if (argc >= 3) DISPLAY = atoi(argv[2]);
  else DISPLAY = 0;

  strcpy(encoder,argv[1]);
  strcpy(decoder,argv[1]);
  strcat(encoder,".enc");
  strcat(decoder,".dec");
  printf("encoder = %s\n",encoder);
  printf("decoder = %s\n",decoder);

  gaussian(fp, &s);
  row_perm(&s);
  fclose(fp);
  if ((fp = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",argv[1]);
    exit(-1);
  }
  read_spmat_file(fp, &t);
  print_result(&s,&t,encoder,decoder);
}
