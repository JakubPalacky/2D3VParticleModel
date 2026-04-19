#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "save.h"

void save_n_E_U(double *a_U_old, INTENZITA *a_E_old, KONCENTRACE **a_conc, KONCENTRACE **a_conc_PO, KONCENTRACE **a_conc_porucha, int a_iterace, int kdy)
{
    int       i,j,k,l,n,v_reduk;
    int       N[TYPY_PO];
    char      filename[20];
    double    v_max[TYPY_PO];
    double    v[TYPY_PO][1000];
    double    v_maximum[TYPY_PO] = {2e6,2e3};
    double    print_velocity, print_velocity_direction;



    sprintf(filename, "potencial_%d.txt", kdy);
    if ((potencial_w     =  fopen(filename,"w")) 		== NULL) { printf("Soubor potencial_%d.txt se nepodarilo otevrit. \n",kdy);}
    sprintf(filename, "koncentrace_%d.txt", kdy);
    if ((koncentrace_w   =  fopen(filename,"w")) 		== NULL) { printf("Soubor koncentrace_%d.txt se nepodarilo otevrit. \n",kdy);}
    if(PORUCHA==1){
        sprintf(filename, "koncentrace_PO_%d.txt", kdy);
        if ((koncentrace_PO   =  fopen(filename,"w")) 		== NULL) { printf("Soubor koncentrace__PO_%d.txt se nepodarilo otevrit. \n",kdy);}
        sprintf(filename, "koncentrace_porucha_%d.txt", kdy);
        if ((koncentrace_porucha   =  fopen(filename,"w")) 	== NULL) { printf("Soubor koncentrace_porucha_%d.txt se nepodarilo otevrit. \n",kdy);}
    }
    sprintf(filename, "intenzita_%d.txt", kdy);
    if ((intenzita_w     =  fopen(filename,"w")) 		== NULL) { printf("Soubor intenzita_%d.txt se nepodarilo otevrit. \n",kdy);}
    sprintf(filename, "rozdeleni_PO_%d.txt", kdy);
    if ((rychlost_PO     =  fopen(filename,"w"))        == NULL) { printf("Soubor rozdeleni_PO_%d.txt se nepodarilo otevrit. \n",kdy);}

    sprintf(filename, "rozdelenix_PO_%d.txt", kdy);
    if ((rychlostx_PO     =  fopen(filename,"w"))        == NULL) { printf("Soubor rozdelenix_PO_%d.txt se nepodarilo otevrit. \n",kdy);}
    sprintf(filename, "rozdeleniy_PO_%d.txt", kdy);
    if ((rychlosty_PO     =  fopen(filename,"w"))        == NULL) { printf("Soubor rozdeleniy_PO_%d.txt se nepodarilo otevrit. \n",kdy);}
    sprintf(filename, "rozdeleniz_PO_%d.txt", kdy);
    if ((rychlostz_PO     =  fopen(filename,"w"))        == NULL) { printf("Soubor rozdeleniz_PO_%d.txt se nepodarilo otevrit. \n",kdy);}

    sprintf(filename, "rozdeleni_zdroj_%d.txt", kdy);
    if ((velocity_source     =  fopen(filename,"w"))        == NULL) { printf("Soubor rozdeleni_PO_%d.txt se nepodarilo otevrit. \n",kdy);}
    if(PORUCHA==1){
        sprintf(filename, "rozdeleni_porucha_%d.txt", kdy);
        if ((rychlost_porucha=  fopen(filename,"w"))        == NULL) { printf("Soubor rozdeleni_porucha_%d.txt se nepodarilo otevrit. \n",kdy);}
    }
    if(POCET_SOND==1){
        sprintf(filename, "radial_velocity_%d.txt", kdy);
        if ((radial_velocity=  fopen(filename,"w"))        == NULL) { printf("Soubor radial_velocity_%d.txt se nepodarilo otevrit. \n",kdy);}
        sprintf(filename, "angular_velocity_%d.txt", kdy);
        if ((angular_velocity=  fopen(filename,"w"))        == NULL) { printf("Soubor angular_velocity_%d.txt se nepodarilo otevrit. \n",kdy);}

    }
    /* rozdeleni rychlosti pro PO */
    fprintf(rychlost_PO,"cas= %i s \n", a_iterace);
    for(n=0;n<1000;n++){
        for(i=0;i<TYPY_PO;i++){
            print_velocity = n*(v_maximum[i]/1000);
            print_velocity_direction = n*((2*v_maximum[i])/1000) - v_maximum[i];

            fprintf(rychlost_PO, "%e %e ",  print_velocity, v_distribution[i][n]);
            fprintf(rychlostx_PO, "%e %e ", print_velocity_direction, vx_distribution[i][n]);
            fprintf(rychlosty_PO, "%e %e ", print_velocity_direction, vy_distribution[i][n]);
            fprintf(rychlostz_PO, "%e %e ", print_velocity_direction, vz_distribution[i][n]);
        }
        fprintf(rychlost_PO, "\n");
        fprintf(rychlostx_PO, "\n");
        fprintf(rychlosty_PO, "\n");
        fprintf(rychlostz_PO, "\n");
    }
    fprintf(radial_velocity,"cas= %i s \n", a_iterace);
    fprintf(angular_velocity,"cas= %i s \n", a_iterace);

    for(j=0;j<30;j++){
        for(i=0;i<TYPY_PO;i++){
            for(n=0;n<40;n++){
                fprintf(radial_velocity, "%e ",v_radial[i][n][j]);
                fprintf(angular_velocity, "%e ",v_angular[i][n][j]);
            }
        }
        fprintf(radial_velocity, "\n");
        fprintf(angular_velocity, "\n");
    }

    /* rozdeleni rychlosti pro zdroj */
    for(i=0;i<TYPY_ZDROJ;i++) {
        v_max[i]=0.0;
        N[i]=0;
        for(n=0;n<1000;n++) v[i][n]=0.0;
    }
    fprintf(velocity_source,"cas= %i s \n", a_iterace);
    for(n=0; n<N_ZDROJ; n++){
        for(i=0;i<TYPY_ZDROJ;i++){
            if((p_zdroj_1[n].typ==i)&&(v_max[p_zdroj_1[n].typ]<p_zdroj_1[n].v_cel))        v_max[p_zdroj_1[n].typ]=p_zdroj_1[n].v_cel;
        }
    }
    for(i=0;i<TYPY_ZDROJ;i++)  v_max[i]/=1000;
    for(n=0; n<N_ZDROJ; n++){
        for(i=0;i<TYPY_ZDROJ;i++){
            if(p_zdroj_1[n].typ==i){
                N[i]++;
                v_reduk = (int)(p_zdroj_1[n].v_cel/v_max[i]);
                if( v_reduk<1000 )  v[i][v_reduk]+=1.0;
            }
        }
    }
    for(i=0;i<TYPY_ZDROJ;i++)  fprintf(velocity_source, "v_max_%i %f ",i,v_max[i]);
    fprintf(velocity_source,"\n");
    for(n=0;n<1000;n++){
        for(i=0;i<TYPY_ZDROJ;i++){
            v[i][n]/=N[i];
            fprintf(velocity_source, "%f %f ",n*v_max[i],v[i][n]);
            v[i][n]=0.0;
        }
        fprintf(velocity_source, "\n");
        for(i=0;i<TYPY_ZDROJ;i++) v[i][n]=0.0;
    }
    fprintf(velocity_source,"\n");


    /* rozdeleni rychlosti poruchy */
    if(PORUCHA==1){
        for(i=0;i<TYPY_PO;i++) {
            v_max[i]=0.0;
            N[i]=0;
            for(n=0;n<1000;n++) v[i][n]=0.0;
        }
        fprintf(rychlost_porucha,"cas= %i s \n", a_iterace);
        for(n=0; n<N_PORUCHA; n++){
            for(i=0;i<TYPY_PO;i++){
                if((p_porucha[n].typ==i)&&(v_max[p_porucha[n].typ]<p_porucha[n].v_cel))        v_max[p_porucha[n].typ]=p_porucha[n].v_cel;
            }
        }
        for(i=0;i<TYPY_PO;i++)  v_max[i]/=1000;
        for(n=0; n<N_PORUCHA; n++){
            for(i=0;i<TYPY_PO;i++){
                if(p_porucha[n].typ==i){
                    N[i]++;
                    v_reduk = (int)(p_porucha[n].v_cel/v_max[i]);
                    if( v_reduk<1000 )  v[i][v_reduk]+=1.0;
                }
            }
        }
        for(i=0;i<TYPY_PO;i++)  fprintf(rychlost_porucha, "v_max_%i %f ",i,v_max[i]);
        fprintf(rychlost_porucha,"\n");
        for(n=0;n<1000;n++){
            for(i=0;i<TYPY_PO;i++){
                v[i][n]/=N[i];
                fprintf(rychlost_porucha, "%f %f ",n*v_max[i],v[i][n]);
                v[i][n]=0.0;
            }
            fprintf(rychlost_porucha, "\n");
            for(i=0;i<TYPY_PO;i++) v[i][n]=0.0;
        }
        fprintf(rychlost_porucha,"\n");
    }

    /* ukladani dalsich parametru */
    fprintf(koncentrace_w,"cas= %i s \n", a_iterace);
    if(PORUCHA==1){
        fprintf(koncentrace_PO,"cas= %i s \n", a_iterace);
        fprintf(koncentrace_porucha,"cas= %i s \n", a_iterace);
    }
    fprintf(potencial_w,"cas= %i s \n", a_iterace);
    fprintf(intenzita_w,"cas= %i s \n", a_iterace);
    for (k = 0; k < N_ROW; k++) {
        i = k % Nx_UZLY; j = k / Nx_UZLY;
        fprintf(potencial_w, "%d %d %f \n", i, j, (a_U_old[k]));//+(i*h_PO+X0_PO)*EL_POLE) ); // old je jen fiktivni, jde o prvni ulozeni, takze je to aktual
        fprintf(koncentrace_w, "%d %d ", i, j);
        if(PORUCHA==1){
            fprintf(koncentrace_PO, "%d %d ", i, j);
            fprintf(koncentrace_porucha, "%d %d ", i, j);
        }
        for (l = 0; l < TYPY_PO; l++) {
            if(a_iterace==0){
                fprintf(koncentrace_w, "%e ", a_conc[l][k].new);
                if(PORUCHA==1){
                    fprintf(koncentrace_PO, "%e ", a_conc_PO[l][k].new);
                    fprintf(koncentrace_porucha, "%e ", a_conc_porucha[l][k].new);
                }
            }
            else{
                fprintf(koncentrace_w, "%e ", a_conc[l][k].stara);
                if(PORUCHA==1){
                    fprintf(koncentrace_PO, "%e ", a_conc_PO[l][k].stara);
                    fprintf(koncentrace_porucha, "%e ", a_conc_porucha[l][k].stara);
                }
            }
        }
        fprintf(koncentrace_w, "\n");
        if(PORUCHA==1){
            fprintf(koncentrace_PO, "\n");
            fprintf(koncentrace_porucha, "\n");
        }
        fprintf(intenzita_w, "%d %d %f %f %f \n", i, j, a_E_old[k].x, a_E_old[k].y, sqrt(pow(a_E_old[k].x,2)+pow(a_E_old[k].y,2)));
        //fprintf(intenzita_w, "%d %d %f %f %f \n", i, j, a_E_old[k].x, a_E_old[k].y, sqrt(pow(a_E_old[k].x,2)+pow(a_E_old[k].y,2)));
        if ((k>0) && ((k+1)%Nx_UZLY==0)) {
            fprintf(koncentrace_w, "\n");
            if(PORUCHA==1){
                fprintf(koncentrace_PO, "\n");
                fprintf(koncentrace_porucha, "\n");
            }
            fprintf(intenzita_w, "\n");
            fprintf(potencial_w, "\n");
        }
    }
    fclose(koncentrace_w);
    if(PORUCHA==1){
        fclose(koncentrace_PO);
        fclose(koncentrace_porucha);
    }
    fclose(intenzita_w);
    fclose(potencial_w);
    fclose(rychlost_PO);
    fclose(rychlostx_PO);
    fclose(rychlosty_PO);
    fclose(rychlostz_PO);
    if(PORUCHA==1)  fclose(rychlost_porucha);
    if(POCET_SOND==1){
        fclose(radial_velocity);
        fclose(angular_velocity);
    }
}

void save_param(DDCASTICE *save_oblast, DDCASTICE *save_zdroj_1, DDCASTICE *save_zdroj_2, DDCASTICE *save_zdroj_3, DDCASTICE *save_zdroj_4, int N_castic_PO, int N_castic_zdroj, int aktualni_iterace)
{
    size_t navrat;
    if ((oblast =           fopen("oblast.txt","wb"))        == NULL) { printf("Soubor oblast.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_1 =          fopen("zdroj_1.txt", "wb"))      == NULL) { printf("Soubor zdroj_1.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_2 =          fopen("zdroj_2.txt", "wb"))      == NULL) { printf("Soubor zdroj_2.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_3 =          fopen("zdroj_3.txt", "wb"))      == NULL) { printf("Soubor zdroj_3.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_4 =          fopen("zdroj_4.txt", "wb"))      == NULL) { printf("Soubor zdroj_4.txt se nepodarilo otevrit pro zapis. \n");}
    if ((parametry =        fopen("parametry.txt", "w"))     == NULL) { printf("Soubor paramtery.txt se nepodarilo otevrit pro zapis. \n");}
    /* ulozeni parametru pro pripadne cteni */
    fprintf(parametry,"%i \n", aktualni_iterace);
    fprintf(parametry,"%i \n", N_castic_PO);
    fprintf(parametry,"%i \n", N_castic_zdroj);

    /* ulozeni PO */
    navrat = fwrite(save_oblast, sizeof(DDCASTICE), N_castic_PO, oblast);
    if(navrat!=N_castic_PO) {printf("chyba pocet zapsanych castic neodpovida poctu v PO \n");}
    /* ulozeni zdroj_1 */
    navrat = fwrite(save_zdroj_1, sizeof(DDCASTICE), N_castic_zdroj, zdroj_1);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_1 \n");}
    /* ulozeni zdroje_2 */
    navrat = fwrite(save_zdroj_2, sizeof(DDCASTICE), N_castic_zdroj, zdroj_2);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_2 \n");}
    /* ulozeni zdroje_3 */
    navrat = fwrite(save_zdroj_3, sizeof(DDCASTICE), N_castic_zdroj, zdroj_3);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_3 \n");}
    /* ulozeni zdroje_4 */
    navrat = fwrite(save_zdroj_4, sizeof(DDCASTICE), N_castic_zdroj, zdroj_4);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_4 \n");}

    fclose(parametry);
    fclose(oblast);
    fclose(zdroj_1);
    fclose(zdroj_2);
    fclose(zdroj_3);
    fclose(zdroj_4);
}
void save_param_1(DDCASTICE *save_oblast, DDCASTICE *save_zdroj_1, DDCASTICE *save_zdroj_2, DDCASTICE *save_zdroj_3, DDCASTICE *save_zdroj_4, int N_castic_PO, int N_castic_zdroj, int aktualni_iterace)
{
    size_t navrat;
    if ((oblast_1 =           fopen("oblast_1.txt","wb"))        == NULL) { printf("Soubor oblast_1.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_1_1 =          fopen("zdroj_1_1.txt", "wb"))      == NULL) { printf("Soubor zdroj_1_1.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_2_1 =          fopen("zdroj_2_1.txt", "wb"))      == NULL) { printf("Soubor zdroj_2_1.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_3_1 =          fopen("zdroj_3_1.txt", "wb"))      == NULL) { printf("Soubor zdroj_3_1.txt se nepodarilo otevrit pro zapis. \n");}
    if ((zdroj_4_1 =          fopen("zdroj_4_1.txt", "wb"))      == NULL) { printf("Soubor zdroj_4_1.txt se nepodarilo otevrit pro zapis. \n");}
    if ((parametry_1 =        fopen("parametry_1.txt", "w"))     == NULL) { printf("Soubor paramtery_1.txt se nepodarilo otevrit pro zapis. \n");}
    /* ulozeni parametru pro pripadne cteni */
    fprintf(parametry_1,"%i \n", aktualni_iterace);
    fprintf(parametry_1,"%i \n", N_castic_PO);
    fprintf(parametry_1,"%i \n", N_castic_zdroj);

    /* ulozeni PO */
    navrat = fwrite(save_oblast, sizeof(DDCASTICE), N_castic_PO, oblast_1);
    if(navrat!=N_castic_PO) {printf("chyba pocet zapsanych castic neodpovida poctu v PO \n");}
    /* ulozeni zdroj_1 */
    navrat = fwrite(save_zdroj_1, sizeof(DDCASTICE), N_castic_zdroj, zdroj_1_1);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_1 \n");}
    /* ulozeni zdroje_2 */
    navrat = fwrite(save_zdroj_2, sizeof(DDCASTICE), N_castic_zdroj, zdroj_2_1);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_2 \n");}
    /* ulozeni zdroje_3 */
    navrat = fwrite(save_zdroj_3, sizeof(DDCASTICE), N_castic_zdroj, zdroj_3_1);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_3 \n");}
    /* ulozeni zdroje_4 */
    navrat = fwrite(save_zdroj_4, sizeof(DDCASTICE), N_castic_zdroj, zdroj_4_1);
    if(navrat!=N_castic_zdroj) {printf("chyba pocet zapsanych castic neodpovida poctu ve zdroji_4 \n");}

    fclose(parametry_1);
    fclose(oblast_1);
    fclose(zdroj_1_1);
    fclose(zdroj_2_1);
    fclose(zdroj_3_1);
    fclose(zdroj_4_1);
}
void save_model_parameters(FILE *mp_file){
    int i;
    fprintf(mp_file,"This file has all the information conained in global.c and global.h for current run. \n");
    fprintf(mp_file,"N_ITER         %i\n",N_ITER);
    fprintf(mp_file,"IterBetweenSave	%i\n",IterBetweenSave);
    fprintf(mp_file,"CTENI_ZE_SOUBORU	%i\n",CTENI_ZE_SOUBORU);

    fprintf(mp_file,"n0 			%e\n",n0);
    fprintf(mp_file,"n_NEUT		    %e\n",n_NEUT);
    fprintf(mp_file,"Nx_UZLY		%i\n",Nx_UZLY);
    fprintf(mp_file,"X_PO			%e\n",X_PO);
    fprintf(mp_file,"Y_PO			%e\n",Y_PO);
    fprintf(mp_file,"EL_POLE 		%e\n",EL_POLE);
    fprintf(mp_file,"TYPY_PO 		%i\n",TYPY_PO);
    fprintf(mp_file,"TYPY_ZDROJ	    %i\n",TYPY_ZDROJ);
    fprintf(mp_file,"\n POCET_PO ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%i ",POCET_PO[i]);
    fprintf(mp_file,"\n POCET_ZDROJ ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%i ",POCET_ZDROJ[i]);

    fprintf(mp_file,"\n n_NENARUS ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",n_NENARUS[i]);
    fprintf(mp_file,"\n TEPLOTA_PO ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",TEPLOTA_PO[i]);
    fprintf(mp_file,"\n TEPLOTA_ZDROJE ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",TEPLOTA_ZDROJE[i]);
    fprintf(mp_file,"\n hmotnost ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",hmotnost[i]);
    fprintf(mp_file,"\n DT ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",DT[i]);
    fprintf(mp_file,"\n DT_zdroj ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",DT_zdroj[i]);
    fprintf(mp_file,"\n naboj ");
    for(i=0;i<POCET_PO;i++) fprintf(mp_file,"%e ",naboj[i]);

    fprintf(mp_file,"\n POCET_SOND	    %i\n",POCET_SOND);
    fprintf(mp_file,"\n SONDA ");
    for(i=0;i<POCET_SOND;i++) fprintf(mp_file,"%i ",SONDA[i]);
    fprintf(mp_file,"\n R_S ");
    for(i=0;i<POCET_SOND;i++) fprintf(mp_file,"%e ",R_S[i]);
    fprintf(mp_file,"\n X_S ");
    for(i=0;i<POCET_SOND;i++) fprintf(mp_file,"%e ",X_S[i]);
    fprintf(mp_file,"\n Y_S ");
    for(i=0;i<POCET_SOND;i++) fprintf(mp_file,"%e ",Y_S[i]);
    fprintf(mp_file,"\n U_S ");
    for(i=0;i<POCET_SOND;i++) fprintf(mp_file,"%e ",U_S[i]);
}
