/*
  	FILENAME: mkmat.c
  	AUTHOR: Tadashi Wadayama

	NAME: mkmat

	SYNOPSYS: mkmat j k n s
	j: column weight
	k: row weight
	n: code length
	s: seed of random number generator

	DESCRIPTION:

	The program randomly generates a regular sparce matrix.

	EXAMPLE:
	mkmat 3 6 96 1

	HOW TO MAKE:
	gcc -O2 -o mkmat mkmat.c -lm

	HISTORY:
  	SINCE : Mar. 4, 2002
	Copyright (C) Tadashi Wadayama
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int main(int argc,char **argv)
{
  int i,p,q,r,t;			/* loop counter */
  int seed;			/* seed of random number generator */

  int j;			/* comumn weight */
  int k;			/* row weight */
  int n;			/* code length */
  int m;			/* # of check bits */

  int s;

  int* perm;			/* permutation vector */

  int a,b,tmp;	


  if (argc != 5) {
    printf("usage: mkmat j k n seed\n");
    exit(-1);
  }
  j = atoi(argv[1]);
  k = atoi(argv[2]);
  n = atoi(argv[3]);
  seed = atoi(argv[4]);

  //printf("j = %d\n",j);
  //printf("k = %d\n",k);
  //printf("n = %d\n",n);
  //printf("seed = %d\n",seed);

  if (n % k != 0) {
    fprintf(stderr,"n should be a multiple of k.\n");
    exit(-1);
  }
  s = n / k;

  srand48(seed);

  if ((perm = (int*)malloc(sizeof(int)*n)) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }

  m = s*j;

  printf("%d %d\n",n,m);
  printf("%d %d\n",k,j);
  
  for (i = 0; i <= m-1; i++) {
    printf("%d",k);
    if (i  < m-1) printf(" ");
    else printf("\n");
  }

  for (i = 0; i <= n-1; i++) {
    printf("%d",j);
    if (i  < n-1) printf(" ");
    else printf("\n");
  }

	
  for (t = 0; t <= j-1; t++) {
    
    for (i = 0; i <= n-1; i++) perm[i] = i;

    for (i = 0; i <= n * 100; i++) {
      a = floor(drand48() * n);
      b = floor(drand48() * n);
      tmp = perm[a];
      perm[a] = perm[b];
      perm[b] = tmp;
    }

    for (p = 0; p <= s-1; p++) {
      for (q = 0; q <= k-1; q++) {
	for (r = q+1; r <= k-1; r++) {
	  if (perm[p*k+r] < perm[p*k+q]) {
	    tmp = perm[p*k+q];
	    perm[p*k+q] = perm[p*k+r];
	    perm[p*k+r] = tmp;
	  }

	}
      }
    }


    //    for (i = 0; i <= n-1; i++) printf("%d ",perm[i]);printf("\n");
    
    
    for (p = 0; p <= s-1; p++) {
      for (q = 0; q <= k-1; q++) {
	printf("%d",perm[p*k+q]+1);
	if (q < k-1) printf(" ");
	else printf("\n");
      }
      
    }
  }

}
