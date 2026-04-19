#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <pcg_variants.h>

#include "random_pcg.h"
#include "global.h"
#include "init_PO.h"


int N_PO_init;

void INIT_PO(DDCASTICE *p_init)
{
  int i, j, k;
  int N_aktual_init, N_posuv;
  int SPLNENO, v, v_celk;
  double g[2], theta, fi;
  double r_1[POCET_SOND],r_2[POCET_SOND];
  double r_S;
  double c[POCET_SOND], d[POCET_SOND];

  k=0;
  N_posuv = 0;
  N_PO_init = 0;
  for(j=0;j<POCET_SOND;j++){
    r_1[j]= R_S[j];
    r_2[j]= (X_PO < Y_PO) ? X_PO/2.0 : Y_PO/2.0 ;
  }
  for ( i = 0; i < TYPY_PO; i++) {
    N_aktual_init = 0;
    for(j=0;j<POCET_SOND;j++){
        c[j]= -(n_NENARUS[i] * log(r_1[j])) / log(r_2[j]/r_1[j]);
        d[j] = - n_NENARUS[i] / log(r_1[j]/r_2[j]);
    }
    do {
      do  {
        v_celk=0;
        p_init[N_aktual_init+N_posuv].x[0] = -(X_PO / 2) + random(&rng[0]) * X_PO;
        p_init[N_aktual_init+N_posuv].x[1] = -(Y_PO / 2) + random(&rng[0]) * Y_PO;
        for(j=0;j<POCET_SOND;j++){
            if  (pow(p_init[N_aktual_init+N_posuv].x[0]-X_S[j],2)+pow(p_init[N_aktual_init+N_posuv].x[1]-Y_S[j],2)<=pow(R_S[j],2))  v=0;
            else                                                                                                                    v=1;
            if(SONDA[j]==0) v=1; /* kdybych tam sondu nechal ale jen ji zneaktivnil */
            v_celk+=v;
        }
      } while (v_celk!=POCET_SOND);
      /* vklada castice tak aby nebyly uvnitr aktivnich sond, u neaktivnich nevadi */
      SPLNENO=0;
      for(j=0;j<POCET_SOND;j++){
            if(SONDA[j]==1){
                r_S = sqrt(pow(p_init[N_aktual_init+N_posuv].x[0]-X_S[j],2)+pow(p_init[N_aktual_init+N_posuv].x[1]-Y_S[j],2));
                if ( (n_NENARUS[i] * random(&rng[0])) < (c[j] + d[j] * log(r_S)) ) {
                SPLNENO++;
                }
            }
            else SPLNENO++;
      }
      if (SPLNENO == POCET_SOND) {


        /* generovani nahodnych real cisel v intervalu (0,1) */
        /*
        do {g[0]=random(&rng[0]);} while (g[0]==0.0);
        do {g[1]=random(&rng[0]);} while (g[1]==0.0);

        theta=acos(2*random(&rng[0])-1);
        fi=2*M_PI*random(&rng[0]);

        p_init[N_aktual_init+N_posuv].v_cel = vmax[i]*sqrt(-log(g[0])-log(g[1])*pow(cos(2*M_PI*(random(&rng[0]))), 2));
        p_init[N_aktual_init+N_posuv].v[0] = p_init[N_aktual_init+N_posuv].v_cel*sin(theta)*cos(fi);
        p_init[N_aktual_init+N_posuv].v[1] = p_init[N_aktual_init+N_posuv].v_cel*sin(theta)*sin(fi);
        p_init[N_aktual_init+N_posuv].v[2] = p_init[N_aktual_init+N_posuv].v_cel*cos(theta);
        p_init[N_aktual_init+N_posuv].t_s  = (-1)*TAU_MAX[i]*log(1-random(&rng[0]));
        p_init[N_aktual_init+N_posuv].typ  = i;
        */
        if(k == 4*N_ZDROJ) k = 0;
        if(pAllSourcesSteadyVelocity[k].typ == i){
            p_init[N_aktual_init+N_posuv].v_cel = pAllSourcesSteadyVelocity[k].v_cel;
            p_init[N_aktual_init+N_posuv].v[0] = pAllSourcesSteadyVelocity[k].v[0];
            p_init[N_aktual_init+N_posuv].v[1] = pAllSourcesSteadyVelocity[k].v[1];
            p_init[N_aktual_init+N_posuv].v[2] = pAllSourcesSteadyVelocity[k].v[2];
            p_init[N_aktual_init+N_posuv].t_s  = pAllSourcesSteadyVelocity[k].t_s;
            p_init[N_aktual_init+N_posuv].typ  = i;
            k++;
        }
        else{
            do{
                k++;
                if(k == 4*N_ZDROJ) k=0;
                //printf("%i\n",pAllSourcesSteadyVelocity[k].typ);
            }while(pAllSourcesSteadyVelocity[k].typ != i);
            p_init[N_aktual_init+N_posuv].v_cel = pAllSourcesSteadyVelocity[k].v_cel;
            p_init[N_aktual_init+N_posuv].v[0] = pAllSourcesSteadyVelocity[k].v[0];
            p_init[N_aktual_init+N_posuv].v[1] = pAllSourcesSteadyVelocity[k].v[1];
            p_init[N_aktual_init+N_posuv].v[2] = pAllSourcesSteadyVelocity[k].v[2];
            p_init[N_aktual_init+N_posuv].t_s  = pAllSourcesSteadyVelocity[k].t_s;
            p_init[N_aktual_init+N_posuv].typ  = i;
            k++;
        }
        N_aktual_init++;
      }
      N_PO_init++;
    } while (N_aktual_init < POCET_PO[i]);
    N_posuv += POCET_PO[i];
  }

}
