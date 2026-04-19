#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

#include "global.h"
#include "velocity.h"
#include "velocity_porucha.h"

void velocity_porucha(DDCASTICE *castice, INTENZITA *v_E)
{
    double i_notrunc, j_notrunc, dx_h, dy_h, E_x, E_y;
    int k,i,j,r,s,n;

    #pragma omp parallel default(shared) private(i,j,k,s,r,i_notrunc,j_notrunc,dx_h,dy_h,E_x,E_y) num_threads(NT)
    {
        #pragma omp for
        for(n=0; n<N_PORUCHA; n++){
            if(castice[n].in==1){
                /* Vypocet nove rychlosti */
                i_notrunc = (castice[n].x[0] - X0_PO ) / h_PO;
                j_notrunc = (castice[n].x[1] - Y0_PO ) / h_PO;
                i = (int) floor(i_notrunc);
                j = (int) floor(j_notrunc);
                k = i + Nx_UZLY * j;
                dx_h = i_notrunc - i;
                dy_h = j_notrunc - j;
                if ((i<(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                    E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (dx_h) * (1-dy_h) * v_E[k+1].x + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].x + (dx_h) * (dy_h) * v_E[k+1+Nx_UZLY].x;
                    E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (dx_h) * (1-dy_h) * v_E[k+1].y + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].y + (dx_h) * (dy_h) * v_E[k+1+Nx_UZLY].y;
                }
                else if ((i==(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                    E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].x;
                    E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].y;
                }
                else if ((i<(Nx_UZLY-1)) && (j==(Ny_UZLY-1))) {
                    E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (dx_h) * (1-dy_h) * v_E[k+1].x;
                    E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (dx_h) * (1-dy_h) * v_E[k+1].y;
                }
                else {
                    E_x = (1-dx_h) * (1-dy_h) * v_E[k].x;
                    E_y = (1-dx_h) * (1-dy_h) * v_E[k].y;
                }

                castice[n].v[0] = castice[n].v[0] + 1 * KONSTANTA_POHYBU[castice[n].typ]*(E_x+E_extern.x);
                castice[n].v[1] = castice[n].v[1] + 1 * KONSTANTA_POHYBU[castice[n].typ]*(E_y+E_extern.y);
                castice[n].v_cel = sqrt(pow(castice[n].v[0],2) + pow(castice[n].v[1],2) + pow(castice[n].v[2],2));
            }
            else continue;
        }
    }
}
