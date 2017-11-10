#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.h"

#define N 10000000

int main() {

    double sum;
    int i;
    dsfmt_t dsfmt;

    sum = 0.0;

    dsfmt_init_gen_rand(&dsfmt, 7);

    for (i=0; i<N; i++) {
        sum += dsfmt_genrand_close_open(&dsfmt);
    }
    
    printf("Sum of %i random numbers is %f\n", N, sum);

    char* state = dsfmt_state_to_str(&dsfmt, NULL);
    for (i=0; i<10; i++) {
        printf("%i %f\n", i, dsfmt_genrand_close_open(&dsfmt));
    }

    dsfmt_str_to_state(&dsfmt, state, NULL);
    for (i=0; i<10; i++) {
        printf("%i %f\n", i, dsfmt_genrand_close_open(&dsfmt));
    }

    free(state);

    return 0;

}
