#include <math.h>
#include <time.h>
# ifdef _OPENMP
  #include <omp.h>
#endif

#include "global.h"
#include "srazky.h"
#include "Position_adv.h"

void POSITION_ADV(DDCASTICE *p_adv, INTENZITA *v_E, int N_particles)
{
    int n;
    double t_help, x, y, vx, vy, v;
    double i_notrunc, j_notrunc, dx_h, dy_h, E_x, E_y, dt, vex;
    VEKTOR v1help, v2help;
    int k,i,j,r,s;
    int particle_type;

    #pragma omp parallel default(shared) private(t_help, Nr_thread,i,j,k,s,r,i_notrunc,j_notrunc,particle_type,dx_h,dy_h,E_x,E_y,dt) num_threads(NT)
    {
        #pragma omp for
        for(n=0;n<N_particles;n++)
        {
            particle_type = p_adv[n].typ;
            i_notrunc = (p_adv[n].x[0] - X0_PO ) / h_PO;
            j_notrunc = (p_adv[n].x[1] - Y0_PO ) / h_PO;
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

            t_help = p_adv[n].t_s - DT[particle_type];
            // no scatter
            if(t_help > 0.0)
            {
                p_adv[n].x[0] += DT[particle_type] * p_adv[n].v[0] + 0.5 * KONSTANTA_POHYBU[particle_type] * DT[particle_type] * E_x;
                p_adv[n].x[1] += DT[particle_type] * p_adv[n].v[1] + 0.5 * KONSTANTA_POHYBU[particle_type] * DT[particle_type] * E_y;

                p_adv[n].v[0] += KONSTANTA_POHYBU[particle_type] * E_x;
                p_adv[n].v[1] += KONSTANTA_POHYBU[particle_type] * E_y;
                p_adv[n].v_cel = sqrt(pow(p_adv[n].v[0],2) + pow(p_adv[n].v[1],2) + pow(p_adv[n].v[2],2));

                p_adv[n].t_s -= DT[particle_type];
            }
            //scatter
            else
            {
                Nr_thread = omp_get_thread_num();
                do{
                    p_adv[n].x[0] += p_adv[n].t_s * p_adv[n].v[0] + 0.5 * p_adv[n].t_s*p_adv[n].t_s * Qmdiv[particle_type] * E_x;
                    p_adv[n].x[1] += p_adv[n].t_s * p_adv[n].v[1] + 0.5 * p_adv[n].t_s*p_adv[n].t_s * Qmdiv[particle_type] * E_y;

                    p_adv[n].v[0] += p_adv[n].t_s * Qmdiv[particle_type] *E_x;
                    p_adv[n].v[1] += p_adv[n].t_s * Qmdiv[particle_type] *E_y;
                    p_adv[n].v_cel = sqrt(pow(p_adv[n].v[0],2) + pow(p_adv[n].v[1],2) + pow(p_adv[n].v[2],2));

                    srazkove_procesy(p_adv+n, Nr_thread);
                    //printf("1.p_adv[n].t_s %e\n",p_adv[n].t_s);
                    //printf("1.t_help %e\n",t_help);
                    t_help += p_adv[n].t_s;
                    //printf("2.t_help %e\n",t_help);
                } while (t_help <= 0.0);
                t_help = (-1) * (t_help - p_adv[n].t_s);
                //printf("4.t_help %e\n",t_help);
                p_adv[n].x[0] += t_help * p_adv[n].v[0] + 0.5* t_help*t_help * Qmdiv[particle_type] * E_x; //t_help is positive here!
                p_adv[n].x[1] += t_help * p_adv[n].v[1] + 0.5* t_help*t_help * Qmdiv[particle_type] * E_y;

                p_adv[n].v[0] += t_help * Qmdiv[particle_type] *E_x;
                p_adv[n].v[1] += t_help * Qmdiv[particle_type] *E_y;
                p_adv[n].v_cel = sqrt(pow(p_adv[n].v[0],2) + pow(p_adv[n].v[1],2) + pow(p_adv[n].v[2],2));

                p_adv[n].t_s -= t_help;

                //printf("2.p_adv[n].t_s %e\n",p_adv[n].t_s);
            }
        }
    }
}
