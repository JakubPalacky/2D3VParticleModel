#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

#include "global.h"
#include "velocity.h"

void velocity(DDCASTICE *castice, INTENZITA *v_E)
{
    double i_notrunc, j_notrunc, dx_h, dy_h, E_x, E_y, dt, vex;
    VEKTOR v1help, v2help;
    int k,i,j,r,s,n;
    int particle_type;

    #pragma omp parallel default(shared) private(i,j,k,s,r,i_notrunc,j_notrunc,particle_type,dx_h,dy_h,E_x,E_y,dt) num_threads(NT)
    {
        #pragma omp for
        for(n=0; n<N_PO; n++){
            /* Vypocet nove rychlosti */
            particle_type = castice[n].typ;
            i_notrunc = (castice[n].x[0] - X0_PO ) / h_PO;
            j_notrunc = (castice[n].x[1] - Y0_PO ) / h_PO;
            i = (int) floor(i_notrunc);
            j = (int) floor(j_notrunc);
            k = i + Nx_UZLY * j;
            dx_h = i_notrunc - i;
            dy_h = j_notrunc - j;
            if ((i<(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (dx_h) * (1-dy_h) * v_E[k+1].x + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].x + (dx_h) * (dy_h) * v_E[k+1+Nx_UZLY].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (dx_h) * (1-dy_h) * v_E[k+1].y + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].y + (dx_h) * (dy_h) * v_E[k+1+Nx_UZLY].y;
            }
            else if ((i==(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].y;
            }
            else if ((i<(Nx_UZLY-1)) && (j==(Ny_UZLY-1))) {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (dx_h) * (1-dy_h) * v_E[k+1].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (dx_h) * (1-dy_h) * v_E[k+1].y;
            }
            else {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y;
            }
//            if (n<N_PO-N_out){
//                if(MAGNETIC_FIELD){
//                    dt = 0.5;
//                            /* Rotational matrix
//                            B-x (  1,   0,    0;    0, cos, sin;   0, -sin, cos)
//                            B-y (cos,   0, -sin;    0,   1,   0; sin,    0, cos)
//                            B-z (cos, sin,    0; -sin, cos,   0;   0,    0,   1)
//                            aktualne B-x
//                            */
//                    particle_type = castice[n].typ;
//                    v1help.x = castice[n].v[0] + dt * KONSTANTA_POHYBU[particle_type]*(E_x+E_extern.x);
//                    v1help.y = castice[n].v[1] + dt * KONSTANTA_POHYBU[particle_type]*(E_y+E_extern.y);
//                    v1help.z = castice[n].v[2];
//                    v2help.x = Brot_11[particle_type]*v1help.x+Brot_12[particle_type]*v1help.y+Brot_13[particle_type]*v1help.z;
//                    v2help.y = Brot_21[particle_type]*v1help.x+Brot_22[particle_type]*v1help.y+Brot_23[particle_type]*v1help.z;
//                    v2help.z = Brot_31[particle_type]*v1help.x+Brot_32[particle_type]*v1help.y+Brot_33[particle_type]*v1help.z;
//                    castice[n].v[0] = v2help.x + dt * KONSTANTA_POHYBU[particle_type]*(E_x+E_extern.x);
//                    castice[n].v[1] = v2help.y + dt * KONSTANTA_POHYBU[particle_type]*(E_y+E_extern.y);
//                    castice[n].v[2] = v2help.z;
//                    castice[n].v_cel = sqrt(pow(castice[n].v[0],2) + pow(castice[n].v[1],2) + pow(castice[n].v[2],2));
//                }
//                else{
                    //dt = 1;
                    castice[n].v[0] = castice[n].v[0] + KONSTANTA_POHYBU[particle_type]*(E_x);//+E_extern.x);
                    castice[n].v[1] = castice[n].v[1] + KONSTANTA_POHYBU[particle_type]*(E_y);//+E_extern.y);
                    castice[n].v_cel = sqrt(pow(castice[n].v[0],2) + pow(castice[n].v[1],2) + pow(castice[n].v[2],2));
//                }
//            }
//            else{
//                dt = 0.5;
////                if(MAGNETIC_FIELD) {
////                        /*  Back Rotational matrix
////                        B-x* (1,  0,  0;  0, 1, -1;  0, -1, 1)
////                        B-y* (1,  0, -1;  0, 1,  0; -1,  0, 1)
////                        B-z* (1, -1,  0; -1, 1,  0;  0,  0, 1)
////                        aktualne B-x
////                        */
////                    v2help.x = castice[n].v[0];
////                    v2help.y = castice[n].v[1];
////                    v2help.z = castice[n].v[2];
////                    v1help.x = Brot_11[particle_type]*v2help.x+Brot_12[particle_type]*v2help.y+Brot_13[particle_type]*v2help.z;
////                    v1help.y = Brot_21[particle_type]*v2help.x+Brot_22[particle_type]*v2help.y+(-1)*Brot_23[particle_type]*v2help.z;
////                    v1help.z = Brot_31[particle_type]*v2help.x+(-1)*Brot_32[particle_type]*v2help.y+Brot_33[particle_type]*v2help.z;
////                    castice[n].v[0] = v1help.x + dt * KONSTANTA_POHYBU[particle_type]*(E_x+E_extern.x);
////                    castice[n].v[1] = v1help.y + dt * KONSTANTA_POHYBU[particle_type]*(E_y+E_extern.y);
////                    castice[n].v[2] = v1help.z;
////                    castice[n].v_cel = sqrt(pow(castice[n].v[0],2) + pow(castice[n].v[1],2) + pow(castice[n].v[2],2));
////                }
////                else{
//                    castice[n].v[0] = castice[n].v[0] + dt * KONSTANTA_POHYBU[particle_type]*(E_x+E_extern.x);
//                    castice[n].v[1] = castice[n].v[1] + dt * KONSTANTA_POHYBU[particle_type]*(E_y+E_extern.y);
//                    castice[n].v_cel = sqrt(pow(castice[n].v[0],2) + pow(castice[n].v[1],2) + pow(castice[n].v[2],2));
////                }
//
//            }
        }
    }
}

void velocity_after_read(DDCASTICE *castice, INTENZITA *v_E)
{
    // initialization of velocity to t+1/2 after read .. necessary when different DT before and, they have to be placed in DT_old
    double i_notrunc, j_notrunc, dx_h, dy_h, E_x, E_y, dt;
    double DT_old[TYPY_PO] = {1e-11, 1e-8};
    double Const_of_move_old[TYPY_PO];
    int k,i,j,n;
    int particle_type;

    for (i=0; i<TYPY_PO; i++) {
        Const_of_move_old[i]=naboj[i]*DT_old[i]/hmotnost[i];
    }
    dt = 0.5;

    #pragma omp parallel default(shared) private(i,j,k,i_notrunc,j_notrunc,particle_type,dx_h,dy_h,E_x,E_y) num_threads(NT)
    {
        #pragma omp for
        for(n=0; n<N_PO; n++){
            /* Vypocet nove rychlosti */
            particle_type = castice[n].typ;
            i_notrunc = (castice[n].x[0] - X0_PO ) / h_PO;
            j_notrunc = (castice[n].x[1] - Y0_PO ) / h_PO;
            i = (int) floor(i_notrunc);
            j = (int) floor(j_notrunc);
            k = i + Nx_UZLY * j;
            dx_h = i_notrunc - i;
            dy_h = j_notrunc - j;
            if ((i<(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (dx_h) * (1-dy_h) * v_E[k+1].x + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].x + (dx_h) * (dy_h) * v_E[k+1+Nx_UZLY].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (dx_h) * (1-dy_h) * v_E[k+1].y + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].y + (dx_h) * (dy_h) * v_E[k+1+Nx_UZLY].y;
            }
            else if ((i==(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (1-dx_h) * (dy_h) * v_E[k+Nx_UZLY].y;
            }
            else if ((i<(Nx_UZLY-1)) && (j==(Ny_UZLY-1))) {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x + (dx_h) * (1-dy_h) * v_E[k+1].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y + (dx_h) * (1-dy_h) * v_E[k+1].y;
            }
            else {
                //printf("%i %i %i \n",k,i,j);
                E_x = (1-dx_h) * (1-dy_h) * v_E[k].x;
                E_y = (1-dx_h) * (1-dy_h) * v_E[k].y;
            }

            castice[n].v[0] = castice[n].v[0] - dt * Const_of_move_old[particle_type]*(E_x);//+E_extern.x);
            castice[n].v[1] = castice[n].v[1] - dt * Const_of_move_old[particle_type]*(E_y);//+E_extern.y);

            castice[n].v[0] = castice[n].v[0] + dt * KONSTANTA_POHYBU[particle_type]*(E_x);//+E_extern.x);
            castice[n].v[1] = castice[n].v[1] + dt * KONSTANTA_POHYBU[particle_type]*(E_y);//+E_extern.y);
            castice[n].v_cel = sqrt(pow(castice[n].v[0],2) + pow(castice[n].v[1],2) + pow(castice[n].v[2],2));
        }
    }
}
