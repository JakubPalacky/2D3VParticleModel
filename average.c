#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "global.h"
#include "average.h"

void AVERAGE(double *a_U, double *a_U_old, INTENZITA *a_E,	INTENZITA *a_E_old, KONCENTRACE **a_conc, KONCENTRACE **a_conc_PO, KONCENTRACE **a_conc_porucha, int a_iterace)
{
    int i, j, k, l, n;
    double    v_r_help, v_fi_help,inverted_radius;
    double    v_r[TYPY_PO][40][30],v_fi[TYPY_PO][40][30], v[TYPY_PO][1000], vx[TYPY_PO][1000], vy[TYPY_PO][1000], vz[TYPY_PO][1000];
    double    v_maximum[TYPY_PO] = {2e6,2e3};
    double    v_r_max[TYPY_PO] = {2e6,2e3};
    double    v_r_min[TYPY_PO] = {-2e6,-2e3};
    double    v_fi_min[TYPY_PO]= {-2e6,-2e3};
    double    v_fi_max[TYPY_PO]= {2e6,2e3};
    int       v_reduk,vx_reduk,vy_reduk,vz_reduk, r_reduc, v_r_int, v_fi_int;
    int       N[TYPY_PO];
    int       n_part_r[TYPY_PO][40],n_part_fi[TYPY_PO][40];

    for(k=0;k<TYPY_PO;k++){
        N[k] = 0;
        for(i=0;i<40;i++){
            for(j=0;j<30;j++){
                v_r[k][i][j]   = 0.0;
                v_fi[k][i][j]  = 0.0;
            }
            n_part_r[k][i]   = 0;
            n_part_fi[k][i]   = 0;
        }
    }

    for(i=0;i<TYPY_PO;i++){
        for(j=0;j<1000;j++){
            v[i][j]=0.0;
            vx[i][j]=0.0;
            vy[i][j]=0.0;
            vz[i][j]=0.0;
        }
    }
    for(n=0; n<N_PO; n++){
        for(i=0;i<TYPY_PO;i++){
            if(p_oblast[n].typ==i){
                N[i]++;
                v_reduk = (int)( (p_oblast[n].v_cel/v_maximum[i])*1000 );
                vx_reduk = (int)( ((p_oblast[n].v[0]+v_maximum[i])/(2*v_maximum[i]))*1000 );
                vy_reduk = (int)( ((p_oblast[n].v[1]+v_maximum[i])/(2*v_maximum[i]))*1000 );
                vz_reduk = (int)( ((p_oblast[n].v[2]+v_maximum[i])/(2*v_maximum[i]))*1000 );
                if(v_reduk<0) printf("%e %i\n",p_oblast[n].v_cel,n);
                if( v_reduk<1000 )
                    {v[i][v_reduk]+=1.0;}
                if((vx_reduk>=0)&&(vx_reduk<1000))
                    {vx[i][vx_reduk]+=1.0;}
                if((vy_reduk>=0)&&(vy_reduk<1000))
                    {vy[i][vy_reduk]+=1.0;}
                if((vz_reduk>=0)&&(vz_reduk<1000))
                    {vz[i][vz_reduk]+=1.0;}
                /* radial velocities parameters */
                if(POCET_SOND==1){
                    inverted_radius = 1/sqrt(pow(p_oblast[n].x[0]-X_S[0],2)+pow(p_oblast[n].x[1]-Y_S[0],2));
                    r_reduc = (int) floor((sqrt(pow(p_oblast[n].x[0]-X_S[0],2)+pow(p_oblast[n].x[1]-Y_S[0],2))-R_S[0])/Lambda_D);
                    if(r_reduc<0)printf("r_reduc was less than 0 %i\n", r_reduc);
                    else{
                        if(r_reduc < 40){
                            v_r_help = ((p_oblast[n].x[0]-X_S[0]) * p_oblast[n].v[0] + (p_oblast[n].x[1]-Y_S[0]) * p_oblast[n].v[1])*inverted_radius;
                            v_fi_help = ((p_oblast[n].x[0]-X_S[0]) * p_oblast[n].v[1] + (p_oblast[n].x[1]-Y_S[0]) * p_oblast[n].v[0])*inverted_radius;
                            if((v_r_min[i] <= v_r_help)&&(v_r_max[i] >= v_r_help)){
                                v_r_int = (int) floor(((v_r_help-v_r_min[i])/(v_r_max[i]-v_r_min[i]))*29);
                                v_r[i][r_reduc][v_r_int]+=1.0;
                                n_part_r[i][r_reduc]++;
                            }
                            if((v_fi_min[i] <= v_fi_help)&&(v_fi_max[i] >= v_fi_help)){
                                v_fi_int = (int) floor(((v_fi_help-v_fi_min[i])/(v_fi_max[i]-v_fi_min[i]))*29);
                                v_fi[i][r_reduc][v_fi_int]+=1.0;
                                n_part_fi[i][r_reduc]++;
                            }
                        }
                    }
                }
            }
        }
    }
    for(i=0;i<TYPY_PO;i++){
        for(j=0;j<1000;j++){
            v[i][j]/=N[i];
            vx[i][j]/=N[i];
            vy[i][j]/=N[i];
            vz[i][j]/=N[i];
        }
    }
    for(i=0;i<TYPY_PO;i++){
        for(k=0;k<30;k++){
            for(j=0;j<40;j++){
                if(n_part_r[i][j]!=0)   { v_r[i][j][k]/=n_part_r[i][j]; }
                else                    { v_r[i][j][k] = 0.0; }
                if(n_part_fi[i][j]!=0)  { v_fi[i][j][k]/=n_part_fi[i][j]; }
                else                    { v_fi[i][j][k] = 0.0; }
            }
        }
    }

    if (a_iterace == iterace_ustaleni)
    {
        for (i = 0; i < N_ROW; i++){
            a_U_old[i]=a_U[i];
            a_E_old[i].x=a_E[i].x;
            a_E_old[i].y=a_E[i].y;
            for (j = 0; j < TYPY_PO; j++) {
                a_conc[j][i].stara=a_conc[j][i].new;
                a_conc_PO[j][i].stara=a_conc_PO[j][i].new;
                a_conc_porucha[j][i].stara=a_conc_porucha[j][i].new;
            }
        }
        for(i=0;i<TYPY_PO;i++){
            for(j=0;j<1000;j++){
                v_distribution[i][j]    = v[i][j];
                vx_distribution[i][j]    = vx[i][j];
                vy_distribution[i][j]    = vy[i][j];
                vz_distribution[i][j]    = vz[i][j];
            }
        }
        for(i=0;i<TYPY_PO;i++){
            for(k=0;k<30;k++){
                for(j=0;j<40;j++){
                    v_radial[i][j][k]   = v_r[i][j][k];
                    v_angular[i][j][k]  = v_fi[i][j][k];
                }
            }
        }
        //average_counter=0;
    }
    if ((a_iterace > iterace_ustaleni)&&((a_iterace%KROK_PRUMEROVANI)==0))
    {
        l=(a_iterace-iterace_ustaleni)/KROK_PRUMEROVANI;
        for (i = 0; i < N_ROW; i++)
            {
            a_U_old[i]=(l*a_U_old[i]+a_U[i])/(l+1);
            a_E_old[i].x=(l*a_E_old[i].x+a_E[i].x)/(l+1);
            a_E_old[i].y=(l*a_E_old[i].y+a_E[i].y)/(l+1);
            for (j = 0; j < TYPY_PO; j++) {
                a_conc[j][i].stara=(l*a_conc[j][i].stara+a_conc[j][i].new)/(l+1);
                a_conc_PO[j][i].stara=(l*a_conc_PO[j][i].stara+a_conc_PO[j][i].new)/(l+1);
                a_conc_porucha[j][i].stara=(l*a_conc_porucha[j][i].stara+a_conc_porucha[j][i].new)/(l+1);
            }
        }
        for(i=0;i<TYPY_PO;i++){
            for(j=0;j<1000;j++){
                v_distribution[i][j]    = (l*v_distribution[i][j]+v[i][j])/(l+1);
                vx_distribution[i][j]    = (l*vx_distribution[i][j]+vx[i][j])/(l+1);
                vy_distribution[i][j]    = (l*vy_distribution[i][j]+vy[i][j])/(l+1);
                vz_distribution[i][j]    = (l*vz_distribution[i][j]+vz[i][j])/(l+1);
            }
        }
        for(i=0;i<TYPY_PO;i++){
            for(k=0;k<30;k++){
                for(j=0;j<40;j++){
                    v_radial[i][j][k]   =(l*v_radial[i][j][k]+v_r[i][j][k])/(l+1); ;
                    v_angular[i][j][k]  =(l*v_angular[i][j][k]+v_fi[i][j][k])/(l+1);;
                }
            }
        }
        //average_counter++;
    }
}
