#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pcg_variants.h>

#include "random_pcg.h"
#include "global.h"
#include "srazky.h"

void srazkove_procesy(DDCASTICE *castice, int thread_n)
{

double sigma_e[4], p_e[4], g[2]; /* pole pro nahodna cisla */;
double E, Ek, delta_E, new_E, theta, fi, beta, gama,c, energy_test, y, z;
int i;
int scatter_type;
double v_relative, E_cm_ratio;
double m_reduc, mu, test_value;
double v_help;
DDCASTICE neutral;
double v_cm[3], u[3],v1_ucm[3], v1_cm[3], v2_cm[3], v4_cm[3];
double v1_cm_total, v2_cm_total, v4_cm_total, v1_ucm_total;

    do {g[0]= random(&rng[thread_n]);} while(g[0]==0.0);
    do {g[1]= random(&rng[thread_n]);} while(g[1]==0.0);

    mu = 1 / (hmotnost[castice->typ]+hmotnost[1]);

    theta=acos(2*random(&rng[thread_n])-1);
    fi=2*M_PI*random(&rng[thread_n]);
    neutral.v_cel = sqrt(-log(g[0])-log(g[1])*pow(cos(2*M_PI*random(&rng[thread_n])), 2))*vmax[1];
    neutral.v[0] = neutral.v_cel*sin(theta)*cos(fi);
    neutral.v[1] = neutral.v_cel*sin(theta)*sin(fi);
    neutral.v[2] = neutral.v_cel*cos(theta);

    if ( castice->typ==0 )
    {
        test_value = castice->v_cel;
        v_cm[0] = (castice->v[0]*hmotnost[0]+neutral.v[0]*hmotnost[1])*mu;
        v_cm[1] = (castice->v[1]*hmotnost[0]+neutral.v[1]*hmotnost[1])*mu;
        v_cm[2] = (castice->v[2]*hmotnost[0]+neutral.v[2]*hmotnost[1])*mu;
        //v1_cm[0] = castice->v[0] - v_cm[0];
        //v1_cm[1] = castice->v[1] - v_cm[1];
        //v1_cm[2] = castice->v[2] - v_cm[2];
        //v1_cm_total = sqrt(v1_cm[0]*v1_cm[0]+v1_cm[1]*v1_cm[1]+v1_cm[2]*v1_cm[2]);
        //v2_cm[0] = neutral.v[0] - v_cm[0];
        //v2_cm[1] = neutral.v[1] - v_cm[1];
        //v2_cm[2] = neutral.v[2] - v_cm[2];
        //v2_cm_total = sqrt(v2_cm[0]*v2_cm[0]+v2_cm[1]*v2_cm[1]+v2_cm[2]*v2_cm[2]);

        //E_cm_ratio = (hmotnost[0]*v1_cm_total*v1_cm_total)/(hmotnost[1]*v2_cm_total*v2_cm_total);
        v_relative = sqrt(pow(castice->v[0]-neutral.v[0],2)+pow(castice->v[1]-neutral.v[1],2)+pow(castice->v[2]-neutral.v[2],2));
        //v_relative = sqrt(pow(v1_cm[0]-v2_cm[0],2)+pow(v1_cm[1]-v2_cm[1],2)+pow(v1_cm[2]-v2_cm[2],2));
        m_reduc = (hmotnost[0]*hmotnost[1])/(hmotnost[0]+hmotnost[1]);
        Scattered[0]+=1.0;
        E = 0.5*m_reduc*v_relative*v_relative*JtoEV;

        delta_E=0;
        /* modifikovany (o v_relative a M_max) ucinny prurez pro elektrony */
        sigma_e[0]=cross_section_interpolate(sigma_e_ela,counter_e_ela,E);
        sigma_e[1]=cross_section_interpolate(sigma_e_exc,counter_e_exc,E);
        sigma_e[2]=cross_section_interpolate(sigma_e_ion,counter_e_ion,E);

        sigma_e[0]=v_relative*sigma_e[0]/M_MAX[0];
        sigma_e[1]=v_relative*sigma_e[1]/M_MAX[0];
        sigma_e[2]=v_relative*sigma_e[2]/M_MAX[0];

        g[0]=random(&rng[thread_n]);
        if(g[0]<sigma_e[0]) scatter_type = 1; //elastic
        else {
            if(g[0]<(sigma_e[0]+sigma_e[1])) scatter_type = 2; //excitation 2
            else{
                if(g[0]<(sigma_e[0]+sigma_e[1]+sigma_e[2])) scatter_type = 3; //ionization 3
                else {
                        scatter_type = 4; // zero
                        Scattered_zero+=1.0;
                }
            }
        }
        //printf("scatter type %i\n",scatter_type);
        if(scatter_type!=4){
            if(scatter_type == 1) {
                delta_E=0.0;
                theta = acos((2*random(&rng[thread_n]))-1);
                fi = 2*M_PI*random(&rng[thread_n]);

                v1_cm_total = v_relative*hmotnost[1]/(hmotnost[0]+hmotnost[1]);
                v1_cm[0] = v1_cm_total*sin(theta)*cos(fi);
                v1_cm[1] = v1_cm_total*sin(theta)*sin(fi);
                v1_cm[2] = v1_cm_total*cos(theta);
            }
            else if(scatter_type == 2) {
                delta_E=11.5;
                Ek = (E - delta_E)*EVtoJ;
                v_relative = sqrt(2*Ek/m_reduc);
                v1_cm_total = v_relative*hmotnost[1]/(hmotnost[0]+hmotnost[1]);
                //v1_cm_total = sqrt((2*Ek*hmotnost[1])/(hmotnost[0]*hmotnost[1]+pow(hmotnost[0],2)));
                theta = acos((2*random(&rng[thread_n]))-1);
                fi = 2*M_PI*random(&rng[thread_n]);

                v1_cm[0] = v1_cm_total*sin(theta)*cos(fi);
                v1_cm[1] = v1_cm_total*sin(theta)*sin(fi);
                v1_cm[2] = v1_cm_total*cos(theta);
            }
            else if(scatter_type == 3) {
                delta_E=15.8;
                Ek = (E - delta_E)*random(&rng[thread_n])*EVtoJ; //rand 0 to 1 (the E-delata_E is devided between two elekctrons)

                v_relative = sqrt(2*Ek/m_reduc);
                v1_cm_total = v_relative*hmotnost[1]/(hmotnost[0]+hmotnost[1]);

                //v4_cm_total = sqrt((Ek*(hmotnost[1]-hmotnost[0]))/(hmotnost[0]*hmotnost[1]+pow(hmotnost[0],2)));
                theta = acos((2*random(&rng[thread_n]))-1);
                fi = 2*M_PI*random(&rng[thread_n]);

                v1_cm[0] = v1_cm_total*sin(theta)*cos(fi);
                v1_cm[1] = v1_cm_total*sin(theta)*sin(fi);
                v1_cm[2] = v1_cm_total*cos(theta);
/*
                v4_cm[0] = v4_cm_total*sin(theta)*cos(fi);
                v4_cm[1] = v4_cm_total*sin(theta)*sin(fi);
                v4_cm[2] = v4_cm_total*cos(theta);
                u[0] = v4_cm[0]-v4_cm_total;
                u[1] = v4_cm[1];
                u[2] = v4_cm[2];

                gama = (M_PI/4)*random(&rng[thread_n]);
                beta = 2*M_PI*random(&rng[thread_n]);
                c = sin(gama)*v4_cm_total/(2*cos(gama));
                if(beta<=(M_PI/2)){
                    y = c*cos(beta);
                    z = c*sin(beta);
                }
                else if(beta<=M_PI){
                    beta = beta - (M_PI/2);
                    y = -c*sin(beta);
                    z = c*cos(beta);
                }
                else if(beta<=(3*M_PI/2)){
                    beta = beta - M_PI;
                    y = -c*cos(beta);
                    z = -c*sin(beta);
                }
                else{
                    beta = beta - (3*M_PI/2);
                    y = c*sin(beta);
                    z = -c*cos(beta);
                }

                v1_ucm[0] = v4_cm[0]/2;
                v1_ucm[1] = y;
                v1_ucm[2] = z;

                v1_cm[0] = v1_ucm[0] + u[0];
                v1_cm[1] = v1_ucm[1] + u[1];
                v1_cm[2] = v1_ucm[2] + u[2];
                */
            }

            castice->v[0] = v1_cm[0] + v_cm[0];
            castice->v[1] = v1_cm[1] + v_cm[1];
            castice->v[2] = v1_cm[2] + v_cm[2];
            castice->v_cel = sqrt(pow(castice->v[0],2)+pow(castice->v[1],2)+pow(castice->v[2],2));
            //test_value = castice->v_cel;
        }
        //printf("scatter type: %i\n", scatter_type);
        castice->t_s = (-1)*TAU_MAX[0]*log(1-random(&rng[thread_n]));
        //test_value = castice->t_s;
    }
    else if (castice->typ>=1)                                                   /* jde o iont */
    {
        v_relative = sqrt(pow(castice->v[0]-neutral.v[0],2)+pow(castice->v[1]-neutral.v[1],2)+pow(castice->v[2]-neutral.v[2],2));
        m_reduc = (hmotnost[castice->typ]*hmotnost[1])/(hmotnost[castice->typ]+hmotnost[1]);
        E=(0.5*m_reduc*pow(v_relative,2))*JtoEV;
        Scattered[1]+=1.0;
       // printf("%f %f %f %f po ion \n",castice->v[0],castice->v[1],castice->v[2],castice->v_cel);
        if (random(&rng[thread_n]) > 0.42) /* charge transfer */
        {
            castice->v[0] = neutral.v[0];
            castice->v[1] = neutral.v[1];
            castice->v[2] = neutral.v[2];
            castice->v_cel = neutral.v_cel;
        }
        else /* pruzny rozptyl, 0.45 proto, ze je 45% sance na royptyl a zbytek je charge transfer */
        {
            v_cm[0] = (castice->v[0]+neutral.v[0])/2;
            v_cm[1] = (castice->v[1]+neutral.v[1])/2;
            v_cm[2] = (castice->v[2]+neutral.v[2])/2;

            v1_cm[0] = castice->v[0] - v_cm[0];
            v1_cm[1] = castice->v[1] - v_cm[1];
            v1_cm[2] = castice->v[2] - v_cm[2];

            theta = acos((2*random(&rng[thread_n]))-1);
            fi = 2*M_PI*random(&rng[thread_n]);
            v1_cm_total = sqrt(v1_cm[0]*v1_cm[0]+v1_cm[1]*v1_cm[1]+v1_cm[2]*v1_cm[2]);
            v1_cm[0] = v1_cm_total*sin(theta)*cos(fi);
            v1_cm[1] = v1_cm_total*sin(theta)*sin(fi);
            v1_cm[2] = v1_cm_total*cos(theta);

            castice->v[0] = v1_cm[0] + v_cm[0];
            castice->v[1] = v1_cm[1] + v_cm[1];
            castice->v[2] = v1_cm[2] + v_cm[2];
            castice->v_cel = sqrt(pow(castice->v[0],2)+pow(castice->v[1],2)+pow(castice->v[2],2));
        }
        //printf("%f %f %f %f po ion \n",castice->v[0],castice->v[1],castice->v[2],castice->v_cel);
        castice->t_s=(-1)*TAU_MAX[castice->typ]*log(1-random(&rng[thread_n]));
    }
}

double scatter_electron_majorant()
{
double sigma_e[4], m_maximal[10000]; /* pole pro nahodna cisla */;
double E, output, v_relative, m_reduc;
int i;

    m_reduc = (hmotnost[0]*hmotnost[1])/(hmotnost[0]+hmotnost[1]);
    for(i=0;i<10000;i++){
        m_maximal[i] = 0.0;
    }

    for(i=0;i<10000;i++){
        v_relative = i * 1000;
        E = 0.5*m_reduc*v_relative*v_relative*JtoEV;
        /* modifikovany (o v_relative a M_max) ucinny prurez pro elektrony */
        sigma_e[0]=cross_section_interpolate(sigma_e_ela,counter_e_ela,E);
        sigma_e[1]=cross_section_interpolate(sigma_e_exc,counter_e_exc,E);
        sigma_e[2]=cross_section_interpolate(sigma_e_ion,counter_e_ion,E);
//printf("%e %e %e %e \n",E,sigma_e[0],sigma_e[1],sigma_e[2]);
        m_maximal[i] = v_relative*(sigma_e[0]+sigma_e[1]+sigma_e[2]);
    }
    output = 0;
    for(i=0;i<10000;i++){
        if(output<m_maximal[i]) output = m_maximal[i];
    }
    return output;
}

void open_cross_sections_file(int *counter_el, int *counter_ex, int *counter_ion)
{
    FILE *file;
    bool reading_switch = false;
    int part = 0;
    char string[200], string_parse[29];
    char string_compare[29] = "-----------------------------";

    file = fopen("Scattering_e_Ar.txt","r");
    while(fgets(string,200, file)!=NULL){
        strncpy(string_parse, string, 29);
        if(reading_switch && strcmp(string_parse,string_compare)!=0){
            if(part==1) (*counter_el)+=1;
            if(part==3) (*counter_ex)+=1;
            if(part==5) (*counter_ion)+=1;
        }
        if(strcmp(string_parse,string_compare)==0) {
                reading_switch = !reading_switch;
                part++;
        }
    }
    fclose(file);
}

void read_el_cross_sections(double **sigma_el, double **sigma_ex, double **sigma_ion)
{
    bool reading_switch = false;
    int part = 0;
    FILE *file;

    int row;
    char string[200], string_parse[29];
    char string_compare[29] = "-----------------------------";
    char string_1[11],string_2[12];
    char *chptr, *ptr;

    file = fopen("Scattering_e_Ar.txt","r");
    while(fgets(string,200, file)!=NULL){
        strncpy(string_parse, string, 29);
        if(reading_switch && strcmp(string_parse,string_compare)!=0){
            if(part==1) {
                chptr = string + 1;
                strncpy(string_1, chptr, 11);
                chptr = string + 13;
                strncpy(string_2, chptr, 12);
                sigma_el[row][0] = strtod(string_1,&ptr);
                sigma_el[row][1] = strtod(string_2,&ptr);
            }
            if(part==3) {
                chptr = string + 1;
                strncpy(string_1, chptr, 11);
                chptr = string + 13;
                strncpy(string_2, chptr, 12);
                sigma_ex[row][0] = strtod(string_1,&ptr);
                sigma_ex[row][1] = strtod(string_2,&ptr);
            }
            if(part==5) {
                chptr = string + 1;
                strncpy(string_1, chptr, 11);
                chptr = string + 13;
                strncpy(string_2, chptr, 12);
                sigma_ion[row][0] = strtod(string_1,&ptr);
                sigma_ion[row][1] = strtod(string_2,&ptr);
            }
            row++;
        }
        if(strcmp(string_parse,string_compare)==0) {
                reading_switch = !reading_switch;
                part++;
                row = 0;
        }
    }
    fclose(file);
}

double cross_section_interpolate(double **sigma_array, int array_size, double particle_energy)
{
    double sigma_value, a, b, c1, c2;
    int i;

    //check for negative energy and print exception
    if(particle_energy < 0){printf("negative energy in interpolate!!\n"); exit(3);}
    /* energy lesser or equal than Emax*/
    if(particle_energy <= sigma_array[array_size-1][0]){
        /* energy higher or equal than Emin */
        if(particle_energy >= sigma_array[0][0]){
            for(i=0;i<array_size;i++){
                if(particle_energy==sigma_array[i][0]) {
                    sigma_value=sigma_array[i][1];
                    break;
                }
                if((particle_energy>sigma_array[i][0])&&(particle_energy<sigma_array[i+1][0])){
                    c1 = sigma_array[i+1][1]-sigma_array[i][1];
                    c2 = sigma_array[i+1][0]-sigma_array[i][0];
                    a = c1/c2;
                    b = sigma_array[i][1] - sigma_array[i][0] * a;
                    sigma_value = particle_energy * a + b;
                    break;
                }
            }
        }
        /* else it is zero */
        else sigma_value = 0.0;
    }
    else {
        c1 = sigma_array[array_size-1][1]-sigma_array[array_size-2][1];
        c2 = sigma_array[array_size-1][0]-sigma_array[array_size-2][0];
        a = c1/c2;
        b = sigma_array[array_size-2][1] - sigma_array[array_size-2][0] * a;
        sigma_value = particle_energy * a + b;
        /* check for sigma lesser than 0 */
        if(sigma_value<0) sigma_value=0;
    }
    return sigma_value;
}
