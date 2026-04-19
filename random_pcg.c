#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pcg_variants.h>
#include "random_pcg.h"
#include "global.h"

double random(pcg64_random_t *rng_number)
{
    return (double) pcg64_random_r(rng_number)*Rand_modif;
}
