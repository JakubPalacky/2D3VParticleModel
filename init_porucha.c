#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <pcg_variants.h>

#include "random_pcg.h"
#include "global.h"
#include "init_porucha.h"

void INIT_PORUCHA(DDCASTICE *p_init)
{
 double g[2], theta, fi;
 int i, j;

 /* naplneni maxwellovskym rozdelenim */
 for (j=0; j<TYPY_PO; j++) {
    for (i=0; i<POCET_PORUCHA[j]; i++) {
        do {g[0]=random(&rng[0]);} while(g[0]==0.0);    /* generovani nahodnych real cisel v intervalu (0,1> */
        do {g[1]=random(&rng[0]);} while(g[1]==0.0);

        p_init[i+(j*POCET_PORUCHA[j])].v_cel=sqrt(-log(g[0])-log(g[1])*pow(cos(2*M_PI*(random(&rng[0]))), 2))*vmax[j]; /* rozdìlení rychlostí, jejich rozhozeni do smeru */
        theta=acos(2*random(&rng[0])-1);
        fi=2*M_PI*random(&rng[0]);
        p_init[i+(j*POCET_PORUCHA[j])].v[0]=fabs(p_init[i+(j*POCET_PORUCHA[j])].v_cel*sin(theta)*cos(fi));
        p_init[i+(j*POCET_PORUCHA[j])].v[1]=    p_init[i+(j*POCET_PORUCHA[j])].v_cel*sin(theta)*sin(fi);
        p_init[i+(j*POCET_PORUCHA[j])].v[2]=    p_init[i+(j*POCET_PORUCHA[j])].v_cel*cos(theta);
        p_init[i+(j*POCET_PORUCHA[j])].x[0]=(X0_PO-X_PORUCHA)+X_PORUCHA*random(&rng[0]);
        p_init[i+(j*POCET_PORUCHA[j])].x[1]=            Y0_PO+Y_PORUCHA*random(&rng[0]);
        p_init[i+(j*POCET_PORUCHA[j])].typ = j;
        p_init[i+(j*POCET_PORUCHA[j])].in  = 0;
        p_init[i+(j*POCET_PORUCHA[j])].t_s = (-1)*TAU_MAX[j]*log(1-random(&rng[0]));
        //printf("%f \n",p_init[i+(j*POCET_PORUCHA[j])].x[0]);
    }
 }
}

