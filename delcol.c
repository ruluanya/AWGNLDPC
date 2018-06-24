/*
  	FILENAME: delcol.c  	AUTHOR: Tadashi Wadayama

	NAME:
	delcol

	SYNOPSYS: delcol spmat_file (DISPLAY=0/1)
	
	spmat_file : a parity check matrix in spmat format

	DESCRIPTION:
	
	The program deletes columns which have
	overlap of ones greater than 1.

	EXAMPLE:
	delcol delcol.test

	BUGS:

	RELATED DOCUMENTS or SOURCE CODES:
	../spmat.form

	HOW TO MAKE:
	gcc -O2 -o delcol delcol.c

	HISTORY:
  	SINCE : Jan. 25, 2000
	row link part is fixed on June 16,2000.

*/

#include <stdio.h>

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

void print_spmatrix(sparce_matrix* s,sparce_matrix* normal_form)
{

  int m,i;
  item* p;
  int deleted;
  int* information_pos;

  if ((information_pos = (int*)malloc(sizeof(int)*(s->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  deleted = 0;
  for (m = 0; m <= s->M-1; m++) {
    if (normal_form->num_ones_in_row[m] == 0) {
      deleted++;
      p = (s->start_row_list[m]).right;
      for (i = 1; i <= s->num_ones_in_row[m]; i++) {
	s->num_ones_in_col[p->n]--;
	if (s->num_ones_in_col[p->n] == 0) {
	  fprintf(stderr,"0-culumn occurs.\n");
	  exit(-1);
	}
	p = p->right;
      }
    }
  }
  printf("%d %d\n",s->N,s->M-deleted);
  printf("%d %d\n",s->biggest_num_ones_row,s->biggest_num_ones_col);
  for (i = 0; i <= s->M-1; i++) {
    if (normal_form->num_ones_in_row[i] != 0)  
      printf("%d ",s->num_ones_in_row[i]);
  }
  printf("\n");
  for (i = 0; i <= s->N-1; i++) {
    printf("%d ",s->num_ones_in_col[i]);
  }
  printf("\n");
  for (m = 0; m <= s->M-1; m++) {
    if (normal_form->num_ones_in_row[m] != 0) {
      p = (s->start_row_list[m]).right;
      for (i = 1; i <= s->num_ones_in_row[m]; i++) {
	printf("%d ",(p->n)+1);
	p = p->right;
      }
      printf("\n");
    }
  }

  for (i = 0; i <= s->N-1; i++) information_pos[i] = 1;
  for (m = 0; m <= s->M-1; m++) {
    p = (normal_form->start_row_list[m]).right;
    if (p != NULL) information_pos[p->n] = 0;
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

int overlap_chk(int* tmp, sparce_matrix* s, int col_num)
{
  int n,m;
  item* p;
  int i;
  int overlap;

  for (n = 0; n <= col_num-1; n++) {
    overlap = 0;
    p = (s->start_col_list[n]).down;
    /* printf("%d-th col: ",n); */
    for (i = 1; i <= s->num_ones_in_col[n]; i++) {
      /* printf("(%d,%d)->",p->m, p->n); */
      if (tmp[p->m] == 1) overlap++;
      p = p->down;
    }
    /* printf("%d-th column overlap = %d\n",n,overlap); */
    if (overlap > 1) return -1; 
  }
  return 0;
}

void del_col(sparce_matrix* res, sparce_matrix* org)
{
  int n,i,j,m,l;
  item *p,*q,*last_item,*new_item;
  int* tmp_column;
  int current_col;
  
  if ((tmp_column = (int*)malloc(sizeof(int)*(org->M))) 
      == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
 
				/* Initialization  */
  
  if ((res->num_ones_in_row = (int*)malloc(sizeof(int)*(org->M))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((res->num_ones_in_col = (int*)malloc(sizeof(int)*(org->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((res->start_col_list = (item*)malloc(sizeof(item)*(org->N))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  if ((res->start_row_list = (item*)malloc(sizeof(item)*(org->M))) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  for (n = 0; n <= org->N-1; n++) {
    res->num_ones_in_col[n] = 0;
  }
  for (m = 0; m <= org->M-1; m++) {
    res->num_ones_in_row[m] = org->num_ones_in_row[m];
  }

  current_col = 0;
  for (n = 0; n <= org->N-1; n++) {
    for (i = 0; i <= org->M-1; i++) tmp_column[i] = 0;
    
    p = (org->start_col_list[n]).down;
    if (DISPLAY == 1) printf("%d-th col: ",n);

    for (j = 1; j <= org->num_ones_in_col[n]; j++) {
      tmp_column[p->m] = 1;
      p = p->down;
    }

    if (DISPLAY == 1) {
      printf("tmp_column = ");
      for (j = 0; j <= org->M-1; j++) printf("%d ",tmp_column[j]);
      printf("\n");
    }

    if (overlap_chk(tmp_column,res,current_col)== -1) {
      if (DISPLAY == 1) printf("*** deletion ***\n"); 
      for (l = 0; l <= org->M-1; l++) 
	res->num_ones_in_row[l] -= tmp_column[l];
      goto retry;
    }
				/* making column list */
    last_item = &(res->start_col_list[current_col]);
    for (m = 0; m <= org->M-1; m++) {
      if (tmp_column[m] == 1) {
				/* making new item */
	if ((new_item = (item*)malloc(sizeof(item))) == NULL) {
	  fprintf(stderr,"Can't allocate memory\n");
	  exit(-1);
	}
	new_item->m = m;
	new_item->n = current_col;
	res->num_ones_in_col[current_col]++;
	new_item->right = NULL;
	new_item->down = NULL;
	last_item->down = new_item;
	last_item = new_item;
      }
    }
    current_col++;
  retry:
  }
  res->N = current_col;
  res->M = org->M;

  for (m = 0; m <= res->M-1; m++) {
    res->num_ones_in_row[m] = 0;
  }

				/* making row links */
  for (n = 0; n <= res->N-1; n++) {
    p = (res->start_col_list[n]).down;
    for (i = 1; i <= res->num_ones_in_col[n]; i++) {
      q = &(res->start_row_list[p->m]);
      res->num_ones_in_row[p->m]++;
				/* The next line is modified on June 16,2000 */
      for (j = 1; j <= res->num_ones_in_row[p->m]-1; j++)  q = q->right;
      q->right = p; 
      p = p->down;
      }
  }


  res->biggest_num_ones_row = 0;
  for (m = 0; m <= res->M-1; m++) {
    if (res->num_ones_in_row[m] >= res->biggest_num_ones_row) 
      res->biggest_num_ones_row = res->num_ones_in_row[m];
  }


  res->biggest_num_ones_col = 0;
  for (n = 0; n <= res->N-1; n++) {
    if (res->num_ones_in_col[n] >= res->biggest_num_ones_col) 
      res->biggest_num_ones_col = res->num_ones_in_col[n];
  }
}

int main(int argc,char **argv)
{
  FILE* fp;
  sparce_matrix s,t;
  int i;

  if (argc < 2) {
    printf("usage: delcol spmat_file (DISPLAY=0/1)\n");
    exit(-1);
  }
  if ((fp = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",argv[1]);
    exit(-1);
  }
  if (argc >= 3) DISPLAY = atoi(argv[2]);
  else DISPLAY = 0;

  read_spmat_file(fp, &t);  
  fclose(fp);
  if ((fp = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"Can't open %s.\n",argv[1]);
    exit(-1);
  }
  
  if (DISPLAY == 1) print_spmatrix_in_spmatform(&t);
  del_col(&s, &t);

  if (DISPLAY == 1) printf("RESULT====================\n");
  print_spmatrix_in_spmatform(&s);
}


