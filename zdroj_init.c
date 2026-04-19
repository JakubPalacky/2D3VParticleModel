#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <pcg_variants.h>

#include "random_pcg.h"
#include "global.h"
#include "zdroj_init.h"

void zdroj_init(DDCASTICE *p_zdroj_init, double X_i, double X_f, double Y_i,  double Y_f)
{
    int i,j;
    int N_count_help;
    double theta, fi;
    double g[2];

    N_count_help = 0;
    for (j=0; j<TYPY_ZDROJ; j++) {
        for (i=0; i<POCET_ZDROJ[j]; i++) {
            do {g[0]=random(&rng[0]);} while(g[0]==0.0);    /* generovani nahodnych real cisel v intervalu (0,1> */
            do {g[1]=random(&rng[0]);} while(g[1]==0.0);

            /* rozdelení rychlostí */
            p_zdroj_init[i+N_count_help].v_cel=sqrt(-log(g[0])-log(g[1])*pow(cos(2*M_PI*random(&rng[0])), 2))*vmax[j];

            /* rozhozeni rychlosti do smeru */
            theta=acos(2*random(&rng[0])-1);
            fi=2*M_PI*random(&rng[0]);

            p_zdroj_init[i+N_count_help].v[0]=p_zdroj_init[i+N_count_help].v_cel*sin(theta)*cos(fi);
            p_zdroj_init[i+N_count_help].v[1]=p_zdroj_init[i+N_count_help].v_cel*sin(theta)*sin(fi);
            p_zdroj_init[i+N_count_help].v[2]=p_zdroj_init[i+N_count_help].v_cel*cos(theta);

            /* polohy castic ve zdroji */
            p_zdroj_init[i+N_count_help].x[0]=X_i+(X_f-X_i)* random(&rng[0]);
            p_zdroj_init[i+N_count_help].x[1]=Y_i+(Y_f-Y_i)* random(&rng[0]);
            p_zdroj_init[i+N_count_help].x[2]=Z_i+(Z_f-Z_i)* random(&rng[0]);
            p_zdroj_init[i+N_count_help].typ=j;

            /* nahodna volna draha */
            p_zdroj_init[i+N_count_help].t_s=(-1)*TAU_MAX[j]*log(1-random(&rng[0]));
            //printf("%e \n", p_zdroj_init[i+(j*POCET_ZDROJ[j])].t_s);
        }
        N_count_help+=POCET_ZDROJ[j];
    }
}
