/** 
 * @file dSFMT_state.c 
 * @brief test functions to save and load the state of a dSFMT
 *
 * @author Andrea C G Mennucci (Scuola Normale Superiore)
 *
 * Copyright (C) 2010 Andrea C G Mennucci
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <string.h>

#include "dSFMT.h"


static dsfmt_t dsfmt_internal_state;

static inline double uniform_rnd(){
  return dsfmt_genrand_open_close(&dsfmt_internal_state);
}

int main(int argc, char * argv[])
{
  int i;
  double v[100000];
  char *prefix=NULL;

  dsfmt_init_gen_rand(&dsfmt_internal_state, 189739);

  /* convert 3 states to 3 strings */
  char *s1=dsfmt_state_to_str(&dsfmt_internal_state,prefix);
  for(i=0;i<10000;i++) v[i]=uniform_rnd();

  char *s2=dsfmt_state_to_str(&dsfmt_internal_state,prefix);
  for(i=10000;i<20000;i++) v[i]=uniform_rnd();

  char *s3=dsfmt_state_to_str(&dsfmt_internal_state,prefix);
  for(i=20000;i<30000;i++) v[i]=uniform_rnd();

#ifdef SAVE_THE_STATE /* this was used once, to create the file */
  {
    char filename[100];
    sprintf(filename,"dSFMT_str_state.%d.out.txt",DSFMT_MEXP);
    FILE *F=fopen(filename,"w");
    fputs(s2,F);
    fclose(F);
  }
#endif

  /* convert back 3 strings to 3 states */
  /* s1 */
  char *err=dsfmt_str_to_state(&dsfmt_internal_state,s1,prefix);
  if(err) {
      printf(err);
      exit(1);
    }
  for(i=0;i<10000;i++) 
    if(v[i]!=uniform_rnd()) {
      printf("%s, mismatch at %d\n",argv[0],i);
      exit(1);
    }

  /* s3 */
  err=dsfmt_str_to_state(&dsfmt_internal_state,s3,prefix);
  if(err) {
      printf(err);
      exit(1);
    }
  for(i=20000;i<30000;i++)
    if(v[i]!=uniform_rnd()) {
      printf("%s, mismatch at %d\n",argv[0],i);
      exit(1);
    }

  /* s2 */
  err=dsfmt_str_to_state(&dsfmt_internal_state,s2,prefix);
  if(err) {
      printf(err);
      exit(1);
    }
  for(i=10000;i<20000;i++)
    if(v[i]!=uniform_rnd()) {
      printf("%s, mismatch at %d\n",argv[0],i);
      exit(1);
    }

  /* s2 again */
  err=dsfmt_str_to_state(&dsfmt_internal_state,s2,prefix);
  if(err) {
      printf(err);
      exit(1);
    }
  for(i=10000;i<20000;i++)
    if(v[i]!=uniform_rnd()) {
      printf("%s, mismatch at %d\n",argv[0],i);
      exit(1);
    }

  /* s2 from file */
  {
    char filename[100];
    sprintf(filename,"dSFMT_str_state.%d.out.txt",DSFMT_MEXP);
    FILE *F=fopen(filename,"r");
    assert(F);
    err=dsfmt_file_to_state(&dsfmt_internal_state,F,prefix); 
    if(err) {
      printf(err);
      exit(1);
    }
    fclose(F);
  }
  for(i=10000;i<20000;i++)
    if(v[i]!=uniform_rnd()) {
      printf("%s, mismatch at %d\n",argv[0],i);
      exit(1);
    }

  /* s3 again */
  err=dsfmt_str_to_state(&dsfmt_internal_state,s3,prefix);
  if(err) {
      printf(err);
      exit(1);
    }
  for(i=20000;i<30000;i++)
    if(v[i]!=uniform_rnd()) {
      printf("%s, mismatch at %d\n",argv[0],i);
      exit(1);
    }

  printf("%s , test OK\n",argv[0]);
  
  return 0;
}



