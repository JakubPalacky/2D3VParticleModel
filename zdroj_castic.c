#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
# ifdef _OPENMP
  #include <omp.h>
#endif
#include <pcg_variants.h>

#include "random_pcg.h"
#include "global.h"
#include "zdroj_castic.h"
#include "srazky.h"

void zdroj_castic(DDCASTICE *zdroj, int cislo_zdroje) /* 1 levy, 2 horni, 3 pravy, 4 dolni */
{
 double c, t_help, E_tot[TYPY_PO];
 VEKTOR v1help, v2help;// v_direction[TYPY_PO];
 int i, j, k;
 int N_counter_help, N_out_ion, N_out_el;

 N_counter_help = 0;
 N_out_ion = 0;
 N_out_el = 0;
// v_direction[0].x = 0.0;
// v_direction[0].y = 0.0;
// v_direction[1].x = 0.0;
// v_direction[1].y = 0.0;

 for (j=0; j<TYPY_PO; j++) {
    E_tot[j]=0;
 }
 for (j=0; j<TYPY_PO; j++) {
    //printf("draha %f \n", zdroj[1].draha);
    /* polohy */
    #pragma omp parallel default(shared) private(t_help,k, Nr_thread) num_threads(NT)
    {
        #pragma omp for
        for (i=0; i<POCET_ZDROJ[j]; i++) {
            k=i+N_counter_help;
            t_help = zdroj[k].t_s - DT[zdroj[k].typ];
            // no scatter
            if(t_help > 0.0)
            {
                zdroj[k].x[0] += DT_zdroj[j] * zdroj[k].v[0] + 0.5 * DT_zdroj[j]*DT_zdroj[j] *a_zdroj[j].x;
                zdroj[k].x[1] += DT_zdroj[j] * zdroj[k].v[1] + 0.5 * DT_zdroj[j]*DT_zdroj[j] *a_zdroj[j].y;

                zdroj[k].v[0] += DT_zdroj[j] * a_zdroj[j].x;
                zdroj[k].v[1] += DT_zdroj[j] * a_zdroj[j].y;
                zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));

                zdroj[k].t_s -= DT_zdroj[j];
            }
            //scatter
            else
            {
                Nr_thread = omp_get_thread_num();
                do{
                    zdroj[k].x[0] += zdroj[k].t_s * zdroj[k].v[0] + 0.5 * zdroj[k].t_s*zdroj[k].t_s * a_zdroj[j].x;
                    zdroj[k].x[1] += zdroj[k].t_s * zdroj[k].v[1] + 0.5 * zdroj[k].t_s*zdroj[k].t_s * a_zdroj[j].y;

                    zdroj[k].v[0] += zdroj[k].t_s * a_zdroj[j].x;
                    zdroj[k].v[1] += zdroj[k].t_s * a_zdroj[j].y;
                    zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));

                    srazkove_procesy(zdroj+k, Nr_thread);
                    t_help += zdroj[k].t_s;
                } while (t_help <= 0.0);
                t_help = (-1) * (t_help - zdroj[k].t_s);
                zdroj[k].x[0] += t_help * zdroj[k].v[0] + 0.5 * t_help * t_help * a_zdroj[j].x; //t_help is positive here!
                zdroj[k].x[1] += t_help * zdroj[k].v[1] + 0.5 * t_help * t_help * a_zdroj[j].y;

                zdroj[k].v[0] += t_help * a_zdroj[j].x;
                zdroj[k].v[1] += t_help * a_zdroj[j].y;
                zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));

                zdroj[k].t_s -= t_help;
            }
        }
    }
    /* hranice zdroje */
    #pragma omp parallel default(shared) private(c,k,Nr_thread) num_threads(NT)
    {
        #pragma omp for
        for(i=0; i<POCET_ZDROJ[j]; i++){
            Nr_thread = omp_get_thread_num();

            k=i+N_counter_help;
            /* vnejsi hranice */
            if ( zdroj[k].x[0]<(X_ZDROJ_1) ) {
                c=trunc((X_ZDROJ_1-zdroj[k].x[0])/(X_ZDROJ_2-X_ZDROJ_1));
                zdroj[k].x[0]+=(X_ZDROJ_2-X_ZDROJ_1)*(c+1);
            }
            else if ( zdroj[k].x[0]>(X_ZDROJ_4) ) {
                c=trunc((zdroj[k].x[0]-X_ZDROJ_4)/(X_ZDROJ_4-X_ZDROJ_3));
                zdroj[k].x[0]-=(X_ZDROJ_4-X_ZDROJ_3)*(c+1);
            }
            if ( zdroj[k].x[1]<(Y_ZDROJ_1) ) {
                c=trunc((Y_ZDROJ_1-zdroj[k].x[1])/(Y_ZDROJ_2-Y_ZDROJ_1));
                zdroj[k].x[1]+=(Y_ZDROJ_2-Y_ZDROJ_1)*(c+1);
            }
            else if ( zdroj[k].x[1]>(Y_ZDROJ_4) ) {
                c=trunc((zdroj[k].x[1]-Y_ZDROJ_4)/(Y_ZDROJ_4-Y_ZDROJ_3));
                zdroj[k].x[1]-=(Y_ZDROJ_4-Y_ZDROJ_3)*(c+1);
            }
            /* hranice nahore, dole */
            if((cislo_zdroje==1)||(cislo_zdroje==3)){
                if ( zdroj[k].x[1]<(Y_ZDROJ_2) ) {
                    c=trunc((Y_ZDROJ_2-zdroj[k].x[1])/(Y_ZDROJ_3-Y_ZDROJ_2));
                    zdroj[k].x[1]+=(Y_ZDROJ_3-Y_ZDROJ_2)*(c+1);
                }
                else if ( zdroj[k].x[1]>(Y_ZDROJ_3) ) {
                    c=trunc((zdroj[k].x[1]-Y_ZDROJ_3)/(Y_ZDROJ_3-Y_ZDROJ_2));
                    zdroj[k].x[1]-=(Y_ZDROJ_3-Y_ZDROJ_2)*(c+1);
                }
            }
            /* hranice vlevo, vpravo */
            if((cislo_zdroje==2)||(cislo_zdroje==4)){
                if ( zdroj[k].x[0]<(X_ZDROJ_2) ) {
                    c=trunc((X_ZDROJ_2-zdroj[k].x[0])/(X_ZDROJ_3-X_ZDROJ_2));
                    zdroj[k].x[0]+=(X_ZDROJ_3-X_ZDROJ_2)*(c+1);
                }
                else if ( zdroj[k].x[0]>(X_ZDROJ_3) ) {
                    c=trunc((zdroj[k].x[0]-X_ZDROJ_3)/(X_ZDROJ_3-X_ZDROJ_2));
                    zdroj[k].x[0]-=(X_ZDROJ_3-X_ZDROJ_2)*(c+1);
                }
            }
            /* Z hranice */
            if (zdroj[k].x[2]<Z_i) {
                c=trunc((Z_i-zdroj[k].x[2])/(Z_f-Z_i));
                zdroj[k].x[2]+=(Z_f-Z_i)*(c+1);
            }
            else if (zdroj[k].x[2]>Z_f) {
                c=trunc((zdroj[k].x[2]-Z_f)/(Z_f-Z_i));
                zdroj[k].x[2]-=(Z_f-Z_i)*(c+1);
            }
        }
    }
    for(i=0; i<POCET_ZDROJ[j]; i++){
        k=i+N_counter_help;
        /* vnitrni hranice */
        if ( (zdroj[k].x[0]>X_ZDROJ_2)&&(zdroj[k].x[0]<X_ZDROJ_3)&&(zdroj[k].x[1]>Y_ZDROJ_2)&&(zdroj[k].x[1]<Y_ZDROJ_3) ) {
            /* kopirovani do pracovni oblasti */
            if(zdroj[k].typ == 1){
                N_out_ion++;
            }
            else if(zdroj[k].typ == 0){
                N_out_el++;
            }
            if(t<N_ITER-1){
                p_oblast[N_PO].x[0] = zdroj[k].x[0];
                p_oblast[N_PO].x[1] = zdroj[k].x[1];
                p_oblast[N_PO].v[0] = zdroj[k].v[0];
                p_oblast[N_PO].v[1] = zdroj[k].v[1];
                p_oblast[N_PO].v[2] = zdroj[k].v[2];
                p_oblast[N_PO].v_cel= zdroj[k].v_cel;
                p_oblast[N_PO].typ  = zdroj[k].typ;
                p_oblast[N_PO].t_s=zdroj[k].t_s;
                switch (cislo_zdroje)
                {
                    case 1:
                        {
                            p_oblast[N_PO].x[1] = Y0_PO+Y_PO*random(&rng[Nr_thread]);
                            break;
                        }
                    case 2:
                        {
                            p_oblast[N_PO].x[0] = X0_PO+X_PO*random(&rng[Nr_thread]);
                            break;
                        }
                    case 3:
                        {
                            p_oblast[N_PO].x[1] = Y0_PO+Y_PO*random(&rng[Nr_thread]);
                            break;
                        }
                    case 4:
                        {
                            p_oblast[N_PO].x[0] = X0_PO+X_PO*random(&rng[Nr_thread]);
                            break;
                        }
                }
                N_PO++;
                N_out++;
                if (N_REZERVA_OBLAST==N_PO) {
                    N_REZERVA_OBLAST=N_PO+N_PRIDAVEK_OBLAST;
                    p_oblast = (DDCASTICE *)    realloc(p_oblast, N_REZERVA_OBLAST*sizeof(DDCASTICE));
                }
            }
            switch (cislo_zdroje){
                case 1:
                    c=trunc((zdroj[k].x[0]-X_ZDROJ_2)/(X_ZDROJ_2-X_ZDROJ_1));
                    zdroj[k].x[0]-=(X_ZDROJ_2-X_ZDROJ_1)*(c+1);
                    break;
                case 2:
                    c=trunc((Y_ZDROJ_3-zdroj[k].x[1])/(Y_ZDROJ_4-Y_ZDROJ_3));
                    zdroj[k].x[1]+=(Y_ZDROJ_4-Y_ZDROJ_3)*(c+1);
                    break;
                case 3:
                    c=trunc((X_ZDROJ_3-zdroj[k].x[0])/(X_ZDROJ_4-X_ZDROJ_3));
                    zdroj[k].x[0]+=(X_ZDROJ_4-X_ZDROJ_3)*(c+1);
                    break;
                case 4:
                    c=trunc((zdroj[k].x[1]-Y_ZDROJ_2)/(Y_ZDROJ_2-Y_ZDROJ_1));
                    zdroj[k].x[1]-=(Y_ZDROJ_2-Y_ZDROJ_1)*(c+1);
                    break;
            }
            if(c>0)printf("c: %f\n", c);
        }
    }

    /* rychlosti */
//    #pragma omp parallel default(shared) private(k) num_threads(NT)
//    {
//        #pragma omp for
//        for (i=0; i<POCET_ZDROJ[j]; i++) {
//           k=i+N_counter_help;
//            if(MAGNETIC_FIELD){
//                v1help.x = zdroj[k].v[0] + 0.5*at_zdroj[j].x;
//                v1help.y = zdroj[k].v[1] + 0.5*at_zdroj[j].y;
//                v2help.x = Brot_11[j]*v1help.x+Brot_12[j]*v1help.y;
//                v2help.y = Brot_21[j]*v1help.x+Brot_22[j]*v1help.y;
//                zdroj[k].v[0] = v2help.x + 0.5*at_zdroj[j].x;
//                zdroj[k].v[1] = v2help.y + 0.5*at_zdroj[j].y;
//                zdroj[k].v_cel = sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));
//            }
//            else {
//                zdroj[k].v[0]+=a_zdroj[j].x;
//                zdroj[k].v[1]+=a_zdroj[j].y;
                // z-slozku ani nepridavam ... neni tu ve 2D mam el pole pouze Ex, Ey
//                zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));
//            }
//            v_direction[j].x += zdroj[k].v[0];
//            v_direction[j].y += zdroj[k].v[1];
//        }
//    }
    for(i=0; i<POCET_ZDROJ[j]; i++){
        k=i+N_counter_help;
        E_tot[j] += 0.5*(zdroj[k].v_cel*zdroj[k].v_cel*hmotnost[j]);
    }
    N_counter_help+=POCET_ZDROJ[j];
 }
 //printf("N_out_ion %i N_e %i N_celk %i Zdroj %i \n", N_out_ion, N_out-N_out_ion, N_out, cislo_zdroje);
//    if(cislo_zdroje!=4){
//        fprintf(zdroj_tok,"%e   %e  ", (double) (N_out_el/(X_PO*Z_PO*DT_zdroj[0])), (double) ((N_out_ion)/(X_PO*Z_PO*DT_zdroj[1])));
//    }
//    else {
//        fprintf(zdroj_tok,"%e   %e  \n", (double) (N_out_el/(X_PO*Z_PO*DT_zdroj[0])), (double) ((N_out_ion)/(X_PO*Z_PO*DT_zdroj[1])));
//    }

//    v_direction[0].x /= POCET_ZDROJ[0];
//    v_direction[0].y /= POCET_ZDROJ[0];
//    v_direction[1].x /= POCET_ZDROJ[1];
//    v_direction[1].y /= POCET_ZDROJ[1];
//    printf("zdroj %i e: vx = %e vy = %e I: vx = %e vy = %e \n",cislo_zdroje,v_direction[0].x,v_direction[0].y,v_direction[1].x,v_direction[1].y);
E_tot[0]/=(1.5*POCET_ZDROJ[0]*kB);
E_tot[1]/=(1.5*POCET_ZDROJ[1]*kB);
//printf("T zdroj: %e, I: %e \n", E_tot[0], E_tot[1]);
if(t%IterBetweenSave==(IterBetweenSave-1))fprintf(zdroj_energy,"T zdroj: %e, I: %e \n", E_tot[0], E_tot[1]);
//printf("%i %i\n", N_out_ion, N_out_el);
}

void zdroj_castic_init_distribution(DDCASTICE *zdroj, int cislo_zdroje)
{
     double c, t_help, E_tot[TYPY_PO], v_help;
 VEKTOR v1help, v2help;// v_direction[TYPY_PO];
 int i, j, k;
 int N_counter_help, Scatter_counter;
 Scatter_counter = 0;
 Scattered_zero=0.0;
 N_counter_help = 0;
// Scatter_counter= 0;
// v_direction[0].x = 0.0;
// v_direction[0].y = 0.0;
// v_direction[1].x = 0.0;
// v_direction[1].y = 0.0;

 for (j=0; j<TYPY_PO; j++) {
    E_tot[j]=0;
 }
 //printf("%e %e \n", zdroj[0].x[0], zdroj[0].x[1]);
 for (j=0; j<TYPY_PO; j++) {
    //printf("draha %f \n", zdroj[1].draha);
    /* polohy */
    #pragma omp parallel default(shared) private(t_help,k, Nr_thread) num_threads(NT)
    {
        #pragma omp for
        for (i=0; i<POCET_ZDROJ[j]; i++) {
            k=i+N_counter_help;
            t_help = zdroj[k].t_s - DT_zdroj[j];
            // no scatter
            if(t_help > 0.0)
            {
                zdroj[k].v[0] += DT_zdroj[j] * a_zdroj[j].x;
                zdroj[k].v[1] += DT_zdroj[j] * a_zdroj[j].y;
                zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));
                zdroj[k].t_s -= DT_zdroj[j];
            }
            //scatter
            else
            {
                Nr_thread = omp_get_thread_num();
                do{
                    zdroj[k].v[0] += zdroj[k].t_s * a_zdroj[j].x;
                    zdroj[k].v[1] += zdroj[k].t_s * a_zdroj[j].y;
                    zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));
                    srazkove_procesy(zdroj+k, Nr_thread);
                    Scatter_counter++;
                    t_help += zdroj[k].t_s;
                } while (t_help <= 0.0);
                t_help = (-1) * (t_help - zdroj[k].t_s);
                zdroj[k].v[0] += t_help * a_zdroj[j].x; //t_help is positive here!
                zdroj[k].v[1] += t_help * a_zdroj[j].y;
                zdroj[k].v_cel=sqrt(pow(zdroj[k].v[0],2) + pow(zdroj[k].v[1],2) + pow(zdroj[k].v[2],2));
                zdroj[k].t_s -= t_help;
            }
        }
    }

    for(i=0; i<POCET_ZDROJ[j]; i++){
        k=i+N_counter_help;
        E_tot[j] += 0.5*(zdroj[k].v_cel*zdroj[k].v_cel*hmotnost[j]);
    }
    E_tot[j] /= (1.5*POCET_ZDROJ[j]*kB);
    //printf("scattered %i %i %f %e\n",j,Scatter_counter, Scattered_zero, E_tot[j]);
    N_counter_help+=POCET_ZDROJ[j];
    Scatter_counter =0;
 }
//E_tot[0]/=(1.5*POCET_ZDROJ[0]*kB);
//E_tot[1]/=(1.5*POCET_ZDROJ[1]*kB);
//printf("T zdroj: %e, I: %e \n", E_tot[0], E_tot[1]);
//    v_direction[0].x /= POCET_ZDROJ[0];
//    v_direction[0].y /= POCET_ZDROJ[0];
//    v_direction[1].x /= POCET_ZDROJ[1];
//    v_direction[1].y /= POCET_ZDROJ[1];
//    printf("zdroj %i e: vx = %e vy = %e I: vx = %e vy = %e \n",cislo_zdroje,v_direction[0].x,v_direction[0].y,v_direction[1].x,v_direction[1].y);
}
