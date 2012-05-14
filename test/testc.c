#include <stdio.h>
#include "dSFMT.h"

#define N 10000000

int main() {

    double sum;
    int i;

    sum = 0.0;

    dsfmt_gv_init_gen_rand(7);

    for (i=0; i<N; i++) {
        sum += dsfmt_gv_genrand_close_open();
    }
    
    printf("Sum of %i random numbers is %f\n", N, sum);

    return 0;

}
