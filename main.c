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
#include "hranice.h"
#include "PIC.h"
#include "srazky.h"
#include "zdroj_init.h"
#include "zdroj_castic.h"
#include "init_PO.h"
#include "init_porucha.h"
#include "umfpack.h"
#include "average.h"
#include "velocity.h"
#include "velocity_porucha.h"
#include "save.h"
#include "dopad_sonda.h"
#include "Position_adv.h"

int main(void)
{
  FILE      *soubor, *vyvoj_proudu, *energie_PO, *test_position, *test_zdroj, *vyvoj_naboje, *model_parameters;
  int 		*p_flag,*p_Ai,*p_Ap;
  int 		*vypadle_ind;
  double 	*p_constpot, *p_constpot_init, *p_Ax, *p_b;

  double    vzdalenost, test_value;

  void 		*Symbolic, *Numeric;
  double 	*null_PIC = (double *) NULL;

  int 		N_srazeno, N_PORUCHA_init;
  int 		N_VYPADLE, N_HRANICE, N_DOPLNENE, N_SRAZKY;
  int       N_ofPartType[TYPY_PO];
  int       N_out_e[4], N_out_i[4];
  int 		i,j,k,l,m,n,s;
  double 	X_ZDROJ_pom, Y_ZDROJ_pom; // pomocne pro zdroj
  int       init_iterace, ustaleni, ustaleni_pomoc;
  double    x_porucha_average;
  double    PO_T_e, PO_T_i;
  int       dny, hodiny, minuty, sekundy;

  char escape_char;
  char escape_control = 'y';

  size_t    navrat_h;

  clock_t   start_p, start_z, start_k, start_s, start_r;
  clock_t   end_p, end_z, end_k, end_s, end_r;
  double    cas_p, cas_z, cas_k, cas_s, cas_r;

    char      testfile[20];
    FILE   *rychlost_test, *test, *flag_w;
    double v_max[TYPY_PO]={ 2.5e6, 3e3};
    double    v[TYPY_PO][1000];
    double    v_help_roz[TYPY_PO][1000];
    double v_dilek[TYPY_PO];
    double T_tot[TYPY_PO];
    int    iterace_help = 0;
    int    testkdy;

        for(k=0;k<TYPY_PO;k++){
            for(i=0;i<1000;i++){
                v[k][i] = 0;
            }
        }
        for(k=0;k<TYPY_PO;k++){
            v_dilek[k] = v_max[k]/1000;
        }



  //double    vxe, vye, vze, vxi, vyi, vzi;

  //vxe=0.0; vye=0.0; vze=0.0; vxi=0.0; vyi=0.0; vzi=0.0;

    clock_t start=clock();
    srand ((unsigned)time( NULL ));

    for(k=0;k<TYPY_PO;k++){
        N_ofPartType[k] = 0;
    }
    N_srazeno=0;
    ustaleni=0;
    ustaleni_pomoc=0;
    average_counter=10001; /* average.c see */

    for(k=0;k<TYPY_PO;k++){
        for(i=0;i<POCET_SOND;i++){
            for(j=0;j<20;j++){
                POMOC_PROUD[k][i][j]=0;
                POMOC_PROUD[k][i][j]=0;
            }
            for(j=0;j<DELENI_SONDA;j++){
                POMOC_TOK[k][i][j]=0;
                TOK_SONDY[k][i][j]=0.0;
            }
        }
    }
    //prepare scattering
    counter_e_ela = 0;
    counter_e_exc = 0;
    counter_e_ion = 0;
        /* determine how much data about every process I have */
    open_cross_sections_file(&counter_e_ela, &counter_e_exc, &counter_e_ion);
        /* allocate memory for sigma data */
    sigma_e_ela = (double **)malloc(counter_e_ela * sizeof(double *));
    for (i=0; i<counter_e_ela; i++) sigma_e_ela[i] = (double *)malloc(2 * sizeof(double));
    sigma_e_exc = (double **)malloc(counter_e_exc * sizeof(double *));
    for (i=0; i<counter_e_exc; i++) sigma_e_exc[i] = (double *)malloc(2 * sizeof(double));
    sigma_e_ion = (double **)malloc(counter_e_ion * sizeof(double *));
    for (i=0; i<counter_e_ion; i++) sigma_e_ion[i] = (double *)malloc(2 * sizeof(double));
        /* read data into arrays*/
    read_el_cross_sections(sigma_e_ela, sigma_e_exc, sigma_e_ion);

    for (i = 0; i < TYPY_PO; i++) { M_MAX[i]*=pow(M_Used_au/M_Ar_au,2/3);}
    //sigma_e_tot=1/(n_NEUT*lambda_MAX[0]);
    test_value = scatter_electron_majorant();
    if(test_value > M_MAX[0] ) M_MAX[0] = test_value;
    for (i = 0; i < TYPY_PO; i++) { TAU_MAX[i] = 1.0 / (n_NEUT * M_MAX[i]); }
    Tot_charge_in_space = 0.0;

# ifdef _OPENMP
    int num_threads;
    # pragma omp parallel private(Nr_thread, num_threads) num_threads(NT)
    {
      Nr_thread = omp_get_thread_num();
      num_threads = omp_get_num_threads();
      printf("Thread %2d. z %2d threadu.\n", Nr_thread, num_threads);
    }
    # else
    printf("Tento program neni kompilovan pomoci OpenMP.\n");
# endif

    printf("Rng initialization \n");
    Rand_modif = pow(2,-64);
    # ifdef _OPENMP
        # pragma omp parallel private(Nr_thread) num_threads(NT)
        {
            Nr_thread = omp_get_thread_num();
            pcg64_srandom_r(&rng[Nr_thread], time(NULL), (intptr_t)&rng[Nr_thread]);
        }
        #else
            Nr_thread = 0;
            pcg64_srandom_r(&rng[Nr_thread], time(NULL), (intptr_t)&rng[Nr_thread]);
    #endif

    cas_p=0; cas_z=0; cas_k=0; cas_s=0; cas_r=0;

    init_iterace=0;
    N_ZDROJ=0; for (i=0; i<TYPY_ZDROJ; i++)  N_ZDROJ+=POCET_ZDROJ[i];    // Aktualni pocet castic ve zdroji

    N_PO=0;    for (i=0; i<TYPY_PO; i++)  N_PO+=POCET_PO[i];             // Aktualni pocet castic v PO

    /* extern el.field with -1 for particle source because of inverted boundary p_U*/
    E_extern.x = (-1)*EL_POLE; // -1 because of
    E_extern.y = 0.0;

        N_REZERVA_OBLAST=N_PO+N_PRIDAVEK_OBLAST;	                // aktualni pocet prvku pole p_oblast
    N_out=0;		    					                    // pocet castic, ktere prechazeji do PO pri zavolani zdroje
    Z_PO = N_PO / (2 * n0 * X_PO * Y_PO);
    X0_PO=-X_PO/ 2.0;                                           // souradnice leveho dolniho rohu PO
    Y0_PO=-Y_PO/ 2.0;
    h_PO = X_PO / (Nx_UZLY - 1);                                // Velikost hrany vypocetni bunky pro PIC
    Cell_volume = h_PO*h_PO*Z_PO; //volume of PIC cell
    N_REZERVA_VYPADLE = 50000;

    /* okraje zdroju cislovany zleva doprava - X, odspodu nahoru - Y, a Z okraj */
    Z_i=0.0;
    Z_f=Z_PO;                                                   // rozmery zdroje, vyska stejna jako "vyska PO

    X_ZDROJ_pom=N_ZDROJ / (2* n0 * Y_PO * Z_PO);                // a X dopocitam at koncentrace castice vsemi smery
    Y_ZDROJ_pom=N_ZDROJ / (2* n0 * X_PO * Z_PO);                // a Y dopocitam at koncentrace castice vsemi smery

    X_ZDROJ_1= -X_PO/2-X_ZDROJ_pom;
    X_ZDROJ_2= -X_PO/2;
    X_ZDROJ_3=  X_PO/2;
    X_ZDROJ_4=  X_PO/2+X_ZDROJ_pom;

    Y_ZDROJ_1= -Y_PO/2-Y_ZDROJ_pom;
    Y_ZDROJ_2= -Y_PO/2;
    Y_ZDROJ_3=  Y_PO/2;
    Y_ZDROJ_4=  Y_PO/2+Y_ZDROJ_pom;
    /* konec znaceni okraju zdroju */

    /* inic parametru poruchy */
    if(PORUCHA){
        Y_PORUCHA=Y_PO;
        Z_PORUCHA=Z_PO;
        N_PORUCHA=2*n_porucha*X_PORUCHA*Y_PORUCHA*Z_PORUCHA;
        N_PORUCHA_init=N_PORUCHA;
        for(i=0;i<TYPY_PO;i++) {
            POCET_PORUCHA[i]= (int) floor(N_PORUCHA/2);
        }
        N_PORUCHA=0;
        for(i=0;i<TYPY_PO;i++) {
            N_PORUCHA+=POCET_PORUCHA[i];
        }
        printf("N_PORUCHA = %i \n", N_PORUCHA);
    }

    Lambda_D = sqrt((eps0*kB/pow(naboj[0],2))/((n_NENARUS[0]/TEPLOTA_PO[0])+(n_NENARUS[1]/TEPLOTA_PO[1])));
    for (i=0; i<TYPY_PO; i++) {
        Qmdiv[i]=naboj[i]/hmotnost[i];
        KONSTANTA_POHYBU[i]=Qmdiv[i]*DT[i];

        a_zdroj[i].x=(naboj[i]*E_extern.x)/hmotnost[i];
        a_zdroj[i].y=(naboj[i]*E_extern.y)/hmotnost[i];

        vmax[i]=sqrt(2*kB*TEPLOTA_PO[i]/hmotnost[i]);
        PartForCellInBulk[i] = n_NENARUS[i]*Cell_volume;
        for (j=0; j<POCET_SOND; j++) PROUD_SONDY[i][j]=0;
        Scattered[i]=0.0;
        Omega_p[i]=sqrt((n_NENARUS[i]*pow(naboj[0],2))/(hmotnost[i]*eps0));
    }
    Scattered_zero=0.0;
//Lambda_D = sqrt((eps0*kB/pow(naboj[0],2))/((n_NENARUS[0]/35000)+(n_NENARUS[1]/TEPLOTA_PO[1])));
    if( h_PO > 3.4*Lambda_D ) {printf("Condition deltaX < 3.4x Lambda_D for PIC not fulfilled");}
    if( DT[0] > 2*Omega_p[0] ) {printf("Condition DT < 2x Omega_P for PIC not fulfilled");}
    if(POCET_SOND == 1){
        if( fabs(X_S[0]-X0_PO)       < 20*Lambda_D ) {printf("There is not enough space between probe and boundaries to compute 20x Lambda_D \n");}
        if( fabs(X0_PO+X_PO-X_S[0])  < 20*Lambda_D ) {printf("There is not enough space between probe and boundaries to compute 20x Lambda_D \n");}
        if( fabs(Y_S[0]-Y0_PO)       < 20*Lambda_D ) {printf("There is not enough space between probe and boundaries to compute 20x Lambda_D \n");}
        if( fabs(Y0_PO+Y_PO-Y_S[0])  < 20*Lambda_D ) {printf("There is not enough space between probe and boundaries to compute 20x Lambda_D \n");}
    }
/* yapis globalnich parametru do souboru */
    if((model_parameters = fopen("model_parameters.txt","w")) == NULL) {printf("Problem s otevrenim souboru model_parameters.txt pro zapis \n");}
    save_model_parameters(model_parameters);
    fclose(model_parameters);

    if(!CTENI_ZE_SOUBORU) {
        printf("Nedochazi k nacteni ze souboru. \n");
        /*
        if ((oblast =           fopen("oblast.txt","wb"))        == NULL) { printf("Soubor oblast.txt se nepodarilo otevrit pro zapis. \n");}
        if ((zdroj_1 =          fopen("zdroj_1.txt", "wb"))      == NULL) { printf("Soubor zdroj_1.txt se nepodarilo otevrit pro zapis. \n");}
        if ((zdroj_2 =          fopen("zdroj_2.txt", "wb"))      == NULL) { printf("Soubor zdroj_2.txt se nepodarilo otevrit pro zapis. \n");}
        if ((zdroj_3 =          fopen("zdroj_3.txt", "wb"))      == NULL) { printf("Soubor zdroj_3.txt se nepodarilo otevrit pro zapis. \n");}
        if ((zdroj_4 =          fopen("zdroj_4.txt", "wb"))      == NULL) { printf("Soubor zdroj_4.txt se nepodarilo otevrit pro zapis. \n");}
        if ((parametry =        fopen("parametry.txt", "w"))     == NULL) { printf("Soubor paramtery.txt se nepodarilo otevrit pro zapis. \n");}
        */
    }
    else {
        printf("Nacitani ze souboru. \n");
        if ((oblast =           fopen("oblast.txt","rb"))        == NULL) { printf("Soubor oblast.txt se nepodarilo otevrit pro cteni. \n");}
        if ((zdroj_1 =          fopen("zdroj_1.txt", "rb"))      == NULL) { printf("Soubor zdroj_1.txt se nepodarilo otevrit pro cteni. \n");}
        if ((zdroj_2 =          fopen("zdroj_2.txt", "rb"))      == NULL) { printf("Soubor zdroj_2.txt se nepodarilo otevrit pro cteni. \n");}
        if ((zdroj_3 =          fopen("zdroj_3.txt", "rb"))      == NULL) { printf("Soubor zdroj_3.txt se nepodarilo otevrit pro cteni. \n");}
        if ((zdroj_4 =          fopen("zdroj_4.txt", "rb"))      == NULL) { printf("Soubor zdroj_4.txt se nepodarilo otevrit pro cteni. \n");}
        if ((parametry =        fopen("parametry.txt", "r"))     == NULL) { printf("Soubor paramtery.txt se nepodarilo otevrit pro cteni. \n");}
        fscanf(parametry, "%i", &init_iterace);
        fscanf(parametry, "%i", &N_PO);
        fscanf(parametry, "%i", &N_ZDROJ);
        printf("pocatecni iterace: %i, N_PO: %i, N_zdroj: %i \n", init_iterace, N_PO, N_ZDROJ);
        if(init_iterace>N_ITER) printf("POCATECNI ITERACE JE VETSI NEZ CILOVA!!! \n");
        fclose(parametry);
        if ((parametry =        fopen("parametry.txt", "w"))    == NULL) { printf("Soubor paramtery.txt se nepodarilo otevrit pro zapis. \n");}
    }

    if(MAGNETIC_FIELD){
        printf("mam mag. pole \n");
        for(i=0; i<TYPY_PO; i++){
            //Omega_c[i]=Qmdiv[i]*BFIELD_Z;
            Omega_c[i]=Qmdiv[i]*BFIELD_X;
        }
        /* Rotational matrix
        B-x (  1,   0,    0;    0, cos, sin;   0, -sin, cos)
        B-y (cos,   0, -sin;    0,   1,   0; sin,    0, cos)
        B-z (cos, sin,    0; -sin, cos,   0;   0,    0,   1)
        */
        for(i=0;i<TYPY_PO;i++){
            Brot_11[i] = 1;
            Brot_12[i] = 0;
            Brot_13[i] = 0;
            Brot_21[i] = 0;
            Brot_22[i] = cos(Omega_c[i]*DT[i]);
            Brot_23[i] = sin(Omega_c[i]*DT[i]);
            Brot_31[i] = 0;
            Brot_32[i] = -sin(Omega_c[i]*DT[i]);
            Brot_33[i] = cos(Omega_c[i]*DT[i]);
        }
    }

    /* otevreni txt souboru */
    if ((soubor          =  fopen("proudy.txt","w"))            == NULL) { printf("Soubor proudy.txt se nepodarilo otevrit. \n");}
//    if ((poloha_1        =  fopen("poloha_1.txt","w"))          == NULL) { printf("soubor poloha_1.txt se nepodarilo otevrit. \n");}
//    if ((poloha_10        =  fopen("poloha_10.txt","w"))          == NULL) { printf("soubor poloha_10.txt se nepodarilo otevrit. \n");}
    if ((koncentrace_smer=  fopen("koncentrace_smer.txt","w"))  == NULL) { printf("Soubor koncentrace_smer.txt se nepodarilo otevrit. \n");}
    if ((vyvoj_proudu    =  fopen("vyvoj_proudu.txt","w"))      == NULL) { printf("Soubor vyvoj_proudu.txt se nepodarilo otevrit. \n");}
    if ((tok_sonda       =  fopen("tok_sonda.txt","w"))         == NULL) { printf("Soubor tok_sonda.txt se nepodarilo otevrit. \n");}
    if ((zdroj_energy       =  fopen("zdroj_energy.txt","w"))         == NULL) { printf("Soubor zdroj_energy.txt se nepodarilo otevrit. \n");}
    if ((zdroj_tok       =  fopen("zdroj_tok.txt","w"))         == NULL) { printf("Soubor zdroj_tok.txt se nepodarilo otevrit. \n");}
    if ((po_tok       =  fopen("po_tok.txt","w"))         == NULL) { printf("Soubor po_tok.txt se nepodarilo otevrit. \n");}

    if ((test_position       =  fopen("position_test.txt","w"))         == NULL) { printf("Soubor position_test.txt se nepodarilo otevrit. \n");}
    if ((test_zdroj       =  fopen("zdroj_test.txt","w"))         == NULL) { printf("Soubor zdroj_test.txt se nepodarilo otevrit. \n");}

    if ((vyvoj_naboje       =  fopen("vyvoj_naboje.txt","w"))         == NULL) { printf("Soubor vyvoj_naboje.txt se nepodarilo otevrit. \n");}
//    if ((energie_PO         =  fopen("energie_PO.txt", "w"))     == NULL) { printf("Soubor energie_PO.txt se nepodarilo otevrit pro zapis. \n");}
    if ((po_energie         = fopen("Po_teplota.txt","w"))      == NULL) {printf("SOubor Po_teplota.txt se nepovedlo otevrit pro zapis \n");}

    /* Alokace dynamickych poli */
    printf("Prvni alokace dynamickych poli ... \n");
    if ((p_oblast      = (DDCASTICE *) malloc(N_REZERVA_OBLAST*sizeof(DDCASTICE))) == NULL )   { printf("Malo pameti pro pole p_oblast. \n");              exit(1); }
    if ((p_zdroj_1     = (DDCASTICE *) malloc(N_ZDROJ*sizeof(DDCASTICE)))          == NULL )   { printf("Malo pameti pro pole p_zdroj_in. \n");            exit(1); }
    if ((p_zdroj_2     = (DDCASTICE *) malloc(N_ZDROJ*sizeof(DDCASTICE)))          == NULL )   { printf("Malo pameti pro pole p_zdroj_in. \n");            exit(1); }
    if ((p_zdroj_3     = (DDCASTICE *) malloc(N_ZDROJ*sizeof(DDCASTICE)))          == NULL )   { printf("Malo pameti pro pole p_zdroj_out. \n");           exit(1); }
    if ((p_zdroj_4     = (DDCASTICE *) malloc(N_ZDROJ*sizeof(DDCASTICE)))          == NULL )   { printf("Malo pameti pro pole p_zdroj_out. \n");           exit(1); }
    if ((p_porucha     = (DDCASTICE *) malloc(N_PORUCHA*sizeof(DDCASTICE)))        == NULL )   { printf("Malo pameti pro pole p_porucha. \n");             exit(1); }

    if ((pAllSourcesSteadyVelocity     = (DDCASTICE *) malloc(4*N_ZDROJ*sizeof(DDCASTICE)))          == NULL )   { printf("Malo pameti pro pole pAllSourcesSteady. \n");           exit(1); }

    /* inicializace poruchy */
    printf("inicializace poruchy... \n");
    if(PORUCHA){
        INIT_PORUCHA(p_porucha);
    }

    /* Inicializace zdroje */
    printf("inicializace zdroje... \n");
    if(!CTENI_ZE_SOUBORU) {
        zdroj_init(p_zdroj_1,X_ZDROJ_1,X_ZDROJ_2,Y_ZDROJ_2,Y_ZDROJ_3);
        zdroj_init(p_zdroj_2,X_ZDROJ_2,X_ZDROJ_3,Y_ZDROJ_3,Y_ZDROJ_4);
        zdroj_init(p_zdroj_3,X_ZDROJ_3,X_ZDROJ_4,Y_ZDROJ_2,Y_ZDROJ_3);
        zdroj_init(p_zdroj_4,X_ZDROJ_2,X_ZDROJ_3,Y_ZDROJ_1,Y_ZDROJ_2);
        /* add 10us run of source to get real distribution function if Efield is here */
        if(EL_POLE>0.0){
            t=0;
            do{
                zdroj_castic_init_distribution(p_zdroj_1,1);    /* bezi zdroj strana x*/
                zdroj_castic_init_distribution(p_zdroj_2,2);
                zdroj_castic_init_distribution(p_zdroj_3,3);
                zdroj_castic_init_distribution(p_zdroj_4,4);
                t++;
                if(t%10000==0) {printf("time %e from %e\n",t*DT[0],Initial_source_time);}
            }while(t*DT[0] < Initial_source_time);
        }
    }
    else {
        navrat_h=fread(p_zdroj_1, sizeof(DDCASTICE), N_ZDROJ, zdroj_1);
        if (navrat_h!=N_ZDROJ) {    printf("nenacetly se vsechny castice do zdroje_x \n");  return 10;  }
        fclose(zdroj_1);
        //if ((zdroj_1 =            fopen("zdroj_1.txt", "wb"))        == NULL) { printf("Soubor zdroj_1.txt se nepodarilo otevrit pro zapis. \n");}

        navrat_h=fread(p_zdroj_2, sizeof(DDCASTICE), N_ZDROJ, zdroj_2);
        if (navrat_h!=N_ZDROJ) {    printf("nenacetly se vsechny castice do zdroje_x \n");  return 10;  }
        fclose(zdroj_2);
        //if ((zdroj_2 =            fopen("zdroj_2.txt", "wb"))        == NULL) { printf("Soubor zdroj_2.txt se nepodarilo otevrit pro zapis. \n");}

        navrat_h=fread(p_zdroj_3, sizeof(DDCASTICE), N_ZDROJ, zdroj_3);
        if (navrat_h!=N_ZDROJ) {    printf("nenacetly se vsechny castice do zdroje_x \n");  return 10;  }
        fclose(zdroj_3);
        //if ((zdroj_3 =            fopen("zdroj_3.txt", "wb"))        == NULL) { printf("Soubor zdroj_3.txt se nepodarilo otevrit pro zapis. \n");}

        navrat_h=fread(p_zdroj_4, sizeof(DDCASTICE), N_ZDROJ, zdroj_4);
        if (navrat_h!=N_ZDROJ) {    printf("nenacetly se vsechny castice do zdroje_x \n");  return 10;  }
        fclose(zdroj_4);
        //if ((zdroj_4 =            fopen("zdroj_4.txt", "wb"))        == NULL) { printf("Soubor zdroj_4.txt se nepodarilo otevrit pro zapis. \n");}
    }

    /* Fill velocities of all particles from sources to array to use them for workspace filling */
if(!CTENI_ZE_SOUBORU) {
    for(i=0; i<N_ZDROJ;i++){
        pAllSourcesSteadyVelocity[i].v[0] = p_zdroj_1[i].v[0];
        pAllSourcesSteadyVelocity[i].v[1] = p_zdroj_1[i].v[1];
        pAllSourcesSteadyVelocity[i].v[2] = p_zdroj_1[i].v[2];
        pAllSourcesSteadyVelocity[i].v_cel= p_zdroj_1[i].v_cel;
        pAllSourcesSteadyVelocity[i].typ  = p_zdroj_1[i].typ;
        pAllSourcesSteadyVelocity[i].t_s  = p_zdroj_1[i].t_s;
        pAllSourcesSteadyVelocity[i+N_ZDROJ].v[0] = p_zdroj_2[i].v[0];
        pAllSourcesSteadyVelocity[i+N_ZDROJ].v[1] = p_zdroj_2[i].v[1];
        pAllSourcesSteadyVelocity[i+N_ZDROJ].v[2] = p_zdroj_2[i].v[2];
        pAllSourcesSteadyVelocity[i+N_ZDROJ].v_cel= p_zdroj_2[i].v_cel;
        pAllSourcesSteadyVelocity[i+N_ZDROJ].typ  = p_zdroj_2[i].typ;
        pAllSourcesSteadyVelocity[i+N_ZDROJ].t_s  = p_zdroj_2[i].t_s;
        pAllSourcesSteadyVelocity[i+2*N_ZDROJ].v[0] = p_zdroj_3[i].v[0];
        pAllSourcesSteadyVelocity[i+2*N_ZDROJ].v[1] = p_zdroj_3[i].v[1];
        pAllSourcesSteadyVelocity[i+2*N_ZDROJ].v[2] = p_zdroj_3[i].v[2];
        pAllSourcesSteadyVelocity[i+2*N_ZDROJ].v_cel= p_zdroj_3[i].v_cel;
        pAllSourcesSteadyVelocity[i+2*N_ZDROJ].typ  = p_zdroj_3[i].typ;
        pAllSourcesSteadyVelocity[i+2*N_ZDROJ].t_s  = p_zdroj_3[i].t_s;
        pAllSourcesSteadyVelocity[i+3*N_ZDROJ].v[0] = p_zdroj_4[i].v[0];
        pAllSourcesSteadyVelocity[i+3*N_ZDROJ].v[1] = p_zdroj_4[i].v[1];
        pAllSourcesSteadyVelocity[i+3*N_ZDROJ].v[2] = p_zdroj_4[i].v[2];
        pAllSourcesSteadyVelocity[i+3*N_ZDROJ].v_cel= p_zdroj_4[i].v_cel;
        pAllSourcesSteadyVelocity[i+3*N_ZDROJ].typ  = p_zdroj_4[i].typ;
        pAllSourcesSteadyVelocity[i+3*N_ZDROJ].t_s  = p_zdroj_4[i].t_s;
    }
}
    /* Inicializace PO */
//
//    if(!CTENI_ZE_SOUBORU) {
//        printf("Inicializace PO ... \n");
//        INIT_PO(p_oblast);
//    }
//    else {
//        navrat_h=fread(p_oblast, sizeof(DDCASTICE),N_PO, oblast);
//        if (navrat_h!=N_PO) {   printf("nenacetly se vsechny castice do PO \n");    return 1;   }
//        fclose(oblast);
//    }
//
//    /* free not needed data*/
//    free((void *) pAllSourcesSteadyVelocity);

    /* Inicializace PIC */
    printf("Inicializace PIC ...\n");

    if((fmod(Y_PO,h_PO))>1.0e-6) printf("Dojde k chybe v uzlech ve smeru Y!! protoze podil Y_PO a h_PO neni cele cislo!! %e \n", fmod(Y_PO,h_PO));
    Ny_UZLY 		= (int) round(1+(Y_PO/h_PO));                                   		            // number of nodes in Y-directionu !!!! be careful !!!!
    N_ROW 		    = Nx_UZLY * Ny_UZLY;														// number of nodes
    N_NOZERO 		= 5 * (Nx_UZLY-2) * (Ny_UZLY-2) + 2 * (Nx_UZLY + Ny_UZLY - 2);    			// Max count of non-zero elements in A matrix for PIC .. for situation with no probes

    if (PERIOD) { Ny_UZLY_PIC = Ny_UZLY - 1; } // here for period Im not certain
    else { Ny_UZLY_PIC = Ny_UZLY; }

    N_ROW_PIC 	= Nx_UZLY * Ny_UZLY_PIC;
    N_NOZERO_PIC 	= 5 * Nx_UZLY * Ny_UZLY_PIC;

    N_BC_NODES = 2*(Nx_UZLY+Ny_UZLY)-4; //No. of nodes at the edge of workspace

    printf("N_ROW = %d, N_ROW_PIC = %d, Nx_UZLY = %d, Ny_UZLY = %d \n", N_ROW, N_ROW_PIC, Nx_UZLY, Ny_UZLY);

    /* Alokace dynamickych poli pro PIC */
    printf("Alokace poli pro PIC ...\n");
    if ((p_flag = (int *)        					malloc(N_ROW_PIC*sizeof(int))) == NULL )          	{ printf("Malo pameti pro pole p_flag. \n");     			exit(1); }
    if ((p_constpot = (double *) 					malloc(N_ROW_PIC*sizeof(double))) == NULL )       	{ printf("Malo pameti pro pole p_constpot. \n"); 			exit(1); }
    if ((p_constpot_init = (double *) 			malloc(N_ROW_PIC*sizeof(double))) == NULL )       	{ printf("Malo pameti pro pole p_constpot_init. \n"); 			exit(1); } //added 20.10
    if ((p_Ax = (double *)       					malloc(N_NOZERO_PIC*sizeof(double))) == NULL )    	{ printf("Malo pameti pro pole p_Ax. \n");       			exit(1); }
    if ((p_Ai = (int *)          					malloc(N_NOZERO_PIC*sizeof(int))) == NULL )       	{ printf("Malo pameti pro pole p_Ai. \n");       			exit(1); }
    if ((p_Ap = (int *)          					malloc((N_ROW_PIC+1)*sizeof(int))) == NULL )      	{ printf("Malo pameti pro pole p_Ap. \n");       			exit(1); }
    if ((p_b = (double *)        					malloc(N_ROW_PIC*sizeof(double))) == NULL )       	{ printf("Malo pameti pro pole p_b. \n");        			exit(1); }
    if ((p_U = (double *)        					malloc(N_ROW*sizeof(double))) == NULL )       		{ printf("Malo pameti pro pole p_U. \n");        			exit(1); }
    if ((p_gradU_BC = (double *)        			malloc(N_ROW*sizeof(double))) == NULL )       		{ printf("Malo pameti pro pole p_U_BC. \n");        			exit(1); } //added 20.10
    if ((p_U_PIC = (double *)    					malloc(N_ROW_PIC*sizeof(double))) == NULL )       	{ printf("Malo pameti pro pole p_U_PIC. \n");    			exit(1); }
    if ((p_U_old = (double *)            			malloc(N_ROW*sizeof(double))) == NULL )     		{ printf("Malo pameti pro pole p_U_old. \n");    			exit(1); }
    if ((p_E = (INTENZITA *)     					malloc(N_ROW*sizeof(INTENZITA))) == NULL )    		{ printf("Malo pameti pro pole p_E. \n");        			exit(1); }
    if ((p_E_old = (INTENZITA *) 					malloc(N_ROW*sizeof(INTENZITA))) == NULL )    		{ printf("Malo pameti pro pole p_E_old. \n");    			exit(1); }
    if ((vypadle_ind = (int *)   					malloc(N_REZERVA_VYPADLE*sizeof(int))) == NULL )	{ printf("Malo pameti pro pole vypadle_ind. \n");			exit(1); }
    if ((p_U_Bound_pointQ = (double *)        	    malloc(N_BC_NODES*sizeof(double))) == NULL )       		{ printf("Malo pameti pro pole p_U_Bound_pointQ. \n");        			exit(1); }

    printf("Alokace poli pro postprocessing ...\n");
    if ( (conc = (KONCENTRACE **)		 			malloc(TYPY_PO*sizeof(KONCENTRACE *)) ) == NULL ) 	{ printf("Malo pameti pro pole conc."); 	     			exit(1); };
    if ( (conc_PO = (KONCENTRACE **)		 		malloc(TYPY_PO*sizeof(KONCENTRACE *)) ) == NULL ) 	{ printf("Malo pameti pro pole conc_PO."); 	     			exit(1); };
    if ( (conc_porucha = (KONCENTRACE **)		 	malloc(TYPY_PO*sizeof(KONCENTRACE *)) ) == NULL ) 	{ printf("Malo pameti pro pole conc_porucha."); 	     			exit(1); };

    for (i = 0; i < TYPY_PO; i++){
        if ( (conc[i]         = (KONCENTRACE *) 	    malloc(N_ROW*sizeof(KONCENTRACE)) )    == NULL )    { printf("Malo pameti pro %d.radek pole conc.", i);    		exit(1); };
        if ( (conc_PO[i]         = (KONCENTRACE *) 	    malloc(N_ROW*sizeof(KONCENTRACE)) )    == NULL )    { printf("Malo pameti pro %d.radek pole conc_PO.", i);    		exit(1); };
        if ( (conc_porucha[i]         = (KONCENTRACE *) malloc(N_ROW*sizeof(KONCENTRACE)) )    == NULL )    { printf("Malo pameti pro %d.radek pole conc_porucha.", i);    		exit(1); };
    }

    //Real_Qpoint_value = Q_PART_OF_POINTCHARGE();
    //printf("charge of point Q %e %f \n", Real_Qpoint_value,Real_Qpoint_value);
    //fprintf(vyvoj_naboje,"Real_Qpoint_value    %e \n",Real_Qpoint_value);
    /* Priprava poli pro PIC */
    //printf("Setting field of Potential at the edge of workspace for point charge at probe position ...\n");
    //POINT_CHARGE_EDGE_POTENTIAL(p_U_Bound_pointQ);
    printf("Priprava poli flag a constpot ...\n");
    FLAG_CONSTPOT(p_flag, p_constpot);
    /*for(i=0;i<N_ROW;i++){
        p_constpot_init[i] = p_constpot[i];
    }*/

    if ((flag_w     =  fopen("flag.txt","w")) 		== NULL) { printf("Soubor flag.txt se nepodarilo otevrit. \n");}
        for (k = 0; k < N_ROW; k++) {
        i = k % Nx_UZLY; j = k / Nx_UZLY;
        fprintf(flag_w, "%d %d %d \n", i, j, p_flag[k]);//+(i*h_PO+X0_PO)*EL_POLE) ); // old je jen fiktivni, jde o prvni ulozeni, takze je to aktual
        if ((k>0) && ((k+1)%Nx_UZLY==0)) {
            fprintf(flag_w, "\n");
        }
    }
    fclose(flag_w);


    printf("Priprava matic p_Ax, p_Ai, p_Ap ...\n");
    PRIPRAV_A(p_Ax, p_Ai, p_Ap, p_flag, p_constpot);

    printf("Realokace matice p_Ax ... \n");
    if ((p_Ax = (double *) realloc(p_Ax,p_Ap[N_ROW_PIC]*sizeof(double))) == NULL ) { printf("Nelze realokovat pole p_Ax. \n"); exit(1); }
    printf("Realokace matice p_Ai ... \n");
    if ((p_Ai = (int *)    realloc(p_Ai,p_Ap[N_ROW_PIC]*sizeof(int)))    == NULL ) { printf("Nelze realokovat pole p_Ai. \n"); exit(1); }

    printf("Vypocet objektu Symbolic a Numeric ...\n");
    (void) umfpack_di_symbolic (N_ROW_PIC ,N_ROW_PIC, p_Ap, p_Ai, p_Ax, &Symbolic, null_PIC, null_PIC);
    (void) umfpack_di_numeric (p_Ap, p_Ai, p_Ax, Symbolic, &Numeric, null_PIC, null_PIC);
    umfpack_di_free_symbolic(&Symbolic);

    printf("pocet castic v oblasti %i \n", N_PO);

        PRIPRAV_B(p_oblast, p_porucha, p_b, p_constpot, conc, conc_PO, conc_porucha);

        printf("force \n");
        start_s=clock();
        (void) umfpack_di_solve (UMFPACK_At, p_Ap, p_Ai, p_Ax, p_U_PIC, p_b, Numeric, null_PIC, null_PIC);
        for (i = 0; i < N_ROW_PIC; i++){
            p_U[i] = p_U_PIC[i];
        }
//        for (i = 0; i < N_ROW; i++) {
//            p_E[i].x = 0.0;
//            p_E[i].y = 0.0;
//        }

        VYPOCET_E(p_U, p_E);
        end_s=clock();
        cas_s+=(double) (end_s-start_s)/CLOCKS_PER_SEC;

    if ((test     =  fopen("potencial_615.txt","w")) 		== NULL) { printf("Soubor potencial_615.txt se nepodarilo otevrit. \n");}
    fprintf(test,"cas= 10 s \n");
    for (k = 0; k < N_ROW; k++) {
        i = k % Nx_UZLY; j = k / Nx_UZLY;
        fprintf(test, "%d %d %f \n", i, j, (p_U[k]));// // old je jen fiktivni, jde o prvni ulozeni, takze je to aktual
        if ((k>0) && ((k+1)%Nx_UZLY==0)) {
            fprintf(test, "\n");
        }
    }
    fclose(test);

    if(!CTENI_ZE_SOUBORU) {
        printf("Inicializace PO ... \n");
        INIT_PO(p_oblast);
    }
    else {
        navrat_h=fread(p_oblast, sizeof(DDCASTICE),N_PO, oblast);
        if (navrat_h!=N_PO) {   printf("nenacetly se vsechny castice do PO \n");    return 1;   }
        fclose(oblast);
    }

    /* free not needed data*/
    free((void *) pAllSourcesSteadyVelocity);

//upravy 20.10.
        PRIPRAV_B(p_oblast, p_porucha, p_b, p_constpot, conc, conc_PO, conc_porucha);
        (void) umfpack_di_solve (UMFPACK_At, p_Ap, p_Ai, p_Ax, p_U_PIC, p_b, Numeric, null_PIC, null_PIC);
        for (i = 0; i < N_ROW_PIC; i++){
            p_U[i] = p_U_PIC[i];
        }

        VYPOCET_E(p_U, p_E);

//those corrections cannot be used for multiple probes!
       // BOUNDARY_GRADIENT(p_U, p_gradU_BC, p_U_Bound_pointQ, &Tot_charge_in_space); //BC boundary correction
       // CONSTPOT_BOUNDARY_CHANGE(p_constpot_init, p_constpot, p_gradU_BC);
//        PRIPRAV_B(p_oblast, p_porucha, p_b, p_constpot, conc, conc_PO, conc_porucha);
//       (void) umfpack_di_solve (UMFPACK_At, p_Ap, p_Ai, p_Ax, p_U_PIC, p_b, Numeric, null_PIC, null_PIC);
//        for (i = 0; i < N_ROW_PIC; i++){
//            p_U[i] = p_U_PIC[i];
//        }
//        VYPOCET_E(p_U, p_E);
//        BOUNDARY_GRADIENT(p_U, p_gradU_BC, p_U_Bound_pointQ, &Tot_charge_in_space); //BC boundary correction
//        CONSTPOT_BOUNDARY_CHANGE(p_constpot_init, p_constpot, p_gradU_BC);

    //if(!CTENI_ZE_SOUBORU) {
        zdroj_castic(p_zdroj_1,1);    /* bezi zdroj */
        zdroj_castic(p_zdroj_2,2);
        zdroj_castic(p_zdroj_3,3);
        zdroj_castic(p_zdroj_4,4);
    //}

    printf("zacina cas \n");
    /* casova cast */
    t=init_iterace;
    do{
        int vypadla, dilek;
//printf("pocet castic v oblasti %i ze zdroje %i\n", N_PO , N_out);
        if (t%1000==0) printf("cyklus %i \n", t);
        //if (t%10000==0) srand ((unsigned)time( NULL ));
        fprintf(vyvoj_naboje,"%e \n",Tot_charge_in_space);
        if(t%200000==0) {
            if(((t/200000)%2)==0)save_param(p_oblast,p_zdroj_1,p_zdroj_2,p_zdroj_3,p_zdroj_4,N_PO,N_ZDROJ,(N_ITER-1));
            else save_param_1(p_oblast,p_zdroj_1,p_zdroj_2,p_zdroj_3,p_zdroj_4,N_PO,N_ZDROJ,(N_ITER-1));
        }

        x_porucha_average=0.0;
        N_out=0;

        for(k=0;k<TYPY_PO;k++){
            N_ofPartType[k] = 0;
        }

//printf("position \n");
        /* cyklus pres polohy */
        start_p=clock();
        POSITION_ADV(p_oblast, p_E, N_PO);

        //fprintf(poloha_1,"%e %e \n", p_oblast[1].x[0],p_oblast[1].x[1]);
        //fprintf(poloha_10,"%e %e \n", p_oblast[10].x[0],p_oblast[10].x[1]);
        /* cyklus poloh castic poruchy */
        /*
        if((PORUCHA)&&(t>=T_PORUCHA)){
            POSITION_ADV(p_porucha, N_PORUCHA);
        }
        */
        /* konec cyklu pres polohy */
        end_p=clock();
        cas_p+= (double) ( end_p - start_p )/CLOCKS_PER_SEC;
//printf("borders \n");
        /* Kontrola hranic castic v PO */
        for(i=0;i<4;i++){
            N_out_e[i]=0;
            N_out_i[i]=0;
        }

        N_VYPADLE=0;
        for(n=0;n<N_PO;n++){
            vypadla = 0;
            vypadla = HRANICNI_PODMINKY(p_oblast+n);
            if (vypadla != 0) {
                //Make here counter for Ne[0,1,2,3], Ni[0,1,2,3] then compare with "out" from source
                /*
                if(vypadla==1){
                    if(p_oblast[n].x[0]<X0_PO) {
                        if(p_oblast[n].typ==0) N_out_e[0]++;
                        else N_out_i[0]++;
                    }
                    if(p_oblast[n].x[0]>(X0_PO+X_PO)) {
                        if(p_oblast[n].typ==0) N_out_e[2]++;
                        else N_out_i[2]++;
                    }
                    if(p_oblast[n].x[1]<Y0_PO) {
                        if(p_oblast[n].typ==0) N_out_e[3]++;
                        else N_out_i[3]++;
                    }
                    if(p_oblast[n].x[1]>(Y0_PO+Y_PO)) {
                        if(p_oblast[n].typ==0) N_out_e[1]++;
                        else N_out_i[1]++;
                    }
                }
                */
                vypadle_ind[N_VYPADLE] = n;
                N_VYPADLE++;
                if (N_VYPADLE == N_REZERVA_VYPADLE) {
                    N_REZERVA_VYPADLE += N_PRIDAVEK_VYPADLE;
                    if ((vypadle_ind = (int *) realloc(vypadle_ind,N_REZERVA_VYPADLE*sizeof(int))) == NULL ) { printf("Nelze realokovat pole vypadle_ind. \n"); exit(1); }
                }
                for(i=0;i<POCET_SOND;i++){
                    if (vypadla==i+2){
                        POMOC_PROUD[p_oblast[n].typ][i][0]++;
                        if (t>=iterace_ustaleni) {
                            dilek=DOPAD_SONDA(p_oblast+n, i);
                            POMOC_TOK[p_oblast[n].typ][i][dilek]++;
                            PROUD_SONDY[p_oblast[n].typ][i]++;
                        } /* castice dopadla na sondu i */
                    }
                }
            }
        }
        //fprintf(po_tok,"%e %e  %e  %e  %e  %e  %e  %e\n",(double) (N_out_e[0]/(Y_PO*Z_PO*DT[0])), (double) (N_out_i[0]/(Y_PO*Z_PO*DT[1])) ,(double) (N_out_e[1]/(X_PO*Z_PO*DT[0])), (double) (N_out_i[1]/(X_PO*Z_PO*DT[1])) ,(double) (N_out_e[2]/(Y_PO*Z_PO*DT[0])), (double) (N_out_i[2]/(Y_PO*Z_PO*DT[1])) ,(double) (N_out_e[3]/(X_PO*Z_PO*DT[0])), (double) (N_out_i[3]/(X_PO*Z_PO*DT[1])));
        //printf("Out of Workspace e I: Left %i %i Up %i %i Right %i %i Down %i %i \n",N_out_e[0],N_out_i[0],N_out_e[1],N_out_i[1],N_out_e[2],N_out_i[2],N_out_e[3],N_out_i[3]);
        /* preusporadani pole p_oblast, vyplneni mezer (vypadlych indexu) casticemi z konce pole,  N_VYPADLE - pocet vypadlych castic */
        for (i=0; i<N_VYPADLE; i++) {
            while ( vypadle_ind[N_VYPADLE-1] == (N_PO-1) ){
                N_PO--;
                N_VYPADLE--;
            }
            j = vypadle_ind[i];
            p_oblast[j] = p_oblast[N_PO-1];
            N_PO--;
        }
        /* Kontrola hranic pro castice poruchy */
        /*
        if((PORUCHA)&&(t>=T_PORUCHA)){
            N_VYPADLE=0;
            for(n=0; n<N_PORUCHA; n++){
                if(p_porucha[n].in==1){
                    vypadla = 0;
                    vypadla = HRANICNI_PODMINKY(p_porucha+n);
                    if (vypadla != 0) {
                        vypadle_ind[N_VYPADLE] = n;
                        N_VYPADLE++;
                        if (N_VYPADLE == N_REZERVA_VYPADLE) {
                            N_REZERVA_VYPADLE += N_PRIDAVEK_VYPADLE;
                            if ((vypadle_ind = (int *) realloc(vypadle_ind,N_REZERVA_VYPADLE*sizeof(int))) == NULL ) { printf("Nelze realokovat pole vypadle_ind. \n"); exit(1); }
                        }
                        for(i=0;i<POCET_SOND;i++){
                            dilek=DOPAD_SONDA(p_oblast+n, i);
                            POMOC_TOK[p_oblast[n].typ][i][dilek]++;
                            if ((vypadla==i+2) && (t>=iterace_ustaleni)) {PROUD_SONDY[p_porucha[n].typ][i]++;}
                        }
                    }
                }
                else continue;
            }
            for (i=0; i<N_VYPADLE; i++) {
                while ( vypadle_ind[N_VYPADLE-1] == (N_PORUCHA-1) ){
                    N_PORUCHA--;
                    N_VYPADLE--;
                }
                j = vypadle_ind[i];
                p_porucha[j] = p_porucha[N_PORUCHA-1];
                N_PORUCHA--;
            }
        }
        */
        //printf("N_Porucha %i \n", N_PORUCHA);
        /* konec kontroly hranic */
        /* hratky s proudy */

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Changed here to 1000 from 10000 !!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((t+1)%IterBetweenSave==0){
                //printf("current \n");
            for(j=0;j<POCET_SOND;j++){
                for(i=0;i<TYPY_PO;i++){
                    fprintf(vyvoj_proudu,"%e    ",(((t+1)/IterBetweenSave)*DT[i]));
                    //printf("%i \n",POMOC_PROUD[i][j][0]);
                    ARITM_PROUD[i][j]=0.0;
                    /*
                    ARITM_ODCH_PROUD[i][j]=0.0;
                    for(k=0;k<20;k++){*/
                    ARITM_PROUD[i][j]=( (double) POMOC_PROUD[i][j][0])*naboj[i]/(DT[i]*IterBetweenSave*(2 * PI * R_S[j] * Z_PO));
                    /*}
                    ARITM_PROUD[i][j]/=20;
                    for(k=0;k<20;k++){
                        ARITM_ODCH_PROUD[i][j]+=fabs((double)POMOC_PROUD[i][j][k]-ARITM_PROUD[i][j]);
                    }
                    ARITM_ODCH_PROUD[i][j]/=20;
                    */
                    fprintf(vyvoj_proudu,"%f ",ARITM_PROUD[i][j]);
                    //printf("%f  ",ARITM_ODCH_PROUD[i][j]/(ARITM_PROUD[0][j]/100));
                    /*if(ustaleni==0){
                        // zde musi byt konkretni zmena pokud pridam vice typu castic
                        if(j==0||j==2){
                            if((ARITM_ODCH_PROUD[0][j]<(6*ARITM_PROUD[0][j]/100))&&(ARITM_ODCH_PROUD[1][j]<(17*ARITM_PROUD[1][j]/100))) ustaleni_pomoc+=2;
                        }
                        else {
                            if((ARITM_ODCH_PROUD[1][j]<(4*ARITM_PROUD[1][j]/100))) ustaleni_pomoc+=2;
                        }

                        //if((ARITM_ODCH_PROUD[i][j]<(6*ARITM_PROUD[0][j]/100))&&(ARITM_ODCH_PROUD[1][j]<(17*ARITM_PROUD[1][j]/100))) ustaleni_pomoc+=2;

                        if(ustaleni_pomoc==(2*POCET_SOND)) {
                            iterace_ustaleni=t+1;
                            N_ITER=iterace_ustaleni+5001;
                            printf("cas %i N_ITER %i\n",t,N_ITER);
                            ustaleni=1;
                        }
                    }
                    for(k=18;k>=0;k--){
                        POMOC_PROUD[i][j][k+1]=POMOC_PROUD[i][j][k];
                    }
                    */
                    POMOC_PROUD[i][j][0]=0;

                }
//                printf("\n");
            }
            //ustaleni_pomoc=0;
            fprintf(vyvoj_proudu,"\n");
        }

        /* konec hratek */

        /* skok potencialu na sonde*/
/*
        if (t==19999){
            for(i=1;i<TYPY_PO;i++){
                DT[i]=DT[0];
                DT_zdroj[i]=DT_zdroj[0];
                KONSTANTA_POHYBU[i]=(naboj[i]*DT[i])/hmotnost[i];
                at_zdroj[i].x=(naboj[i]*DT_zdroj[i]*E_extern.x)/hmotnost[i];
                at_zdroj[i].y=(naboj[i]*DT_zdroj[i]*E_extern.y)/hmotnost[i];
            }
            U_S[1]=7.0;
            FLAG_CONSTPOT(p_flag, p_constpot);
            PRIPRAV_A(p_Ax, p_Ai, p_Ap, p_flag, p_constpot);
        }
*/
//printf("velocity \n");
/*        start_r=clock();

        velocity(p_oblast, p_E);
        //printf("rychlosti 2 \n");

        //if((PORUCHA)&&(t>=T_PORUCHA)) velocity_porucha(p_porucha,p_E);

        //printf("rychlosti 3 \n");
        end_r=clock(); */
//printf("N_PO %i \n", N_PO);

//printf("conc \n");
        start_k=clock();
        //printf("Priprava prave strany B ...\n");
        PRIPRAV_B(p_oblast, p_porucha, p_b, p_constpot, conc, conc_PO, conc_porucha);
        end_k=clock();
        cas_k+= (double) (end_k-start_k)/CLOCKS_PER_SEC;

//printf("force \n");
        start_s=clock();
        (void) umfpack_di_solve (UMFPACK_At, p_Ap, p_Ai, p_Ax, p_U_PIC, p_b, Numeric, null_PIC, null_PIC);
        for (i = 0; i < N_ROW_PIC; i++){
            p_U[i] = p_U_PIC[i];
        }
//        for (i = 0; i < N_ROW; i++) {
//            p_E[i].x = 0.0;
//            p_E[i].y = 0.0;
//        }
        VYPOCET_E(p_U, p_E);
        //BOUNDARY_GRADIENT(p_U, p_gradU_BC,p_U_Bound_pointQ, &Tot_charge_in_space); //BC boundary correction
        //CONSTPOT_BOUNDARY_CHANGE(p_constpot_init, p_constpot, p_gradU_BC);
        end_s=clock();
        cas_s+=(double) (end_s-start_s)/CLOCKS_PER_SEC;
//those corrections cannot be used for multiple probes!

 // test zdroje tady
// printf("source \n");
        start_z=clock();
        zdroj_castic(p_zdroj_1,1);    /* bezi zdroj */
        zdroj_castic(p_zdroj_2,2);
        zdroj_castic(p_zdroj_3,3);
        zdroj_castic(p_zdroj_4,4);
//printf("N_PO %i N_out %i \n", N_PO, N_out);

        end_z=clock();
        cas_z+=(double) (end_z-start_z)/CLOCKS_PER_SEC;

//printf("average \n");
        cas_r+=(double) (end_r-start_r)/CLOCKS_PER_SEC;
        if(t%IterBetweenSave==0) iterace_ustaleni=t; // 100000 or 10000
        AVERAGE(p_U, p_U_old, p_E, p_E_old, conc,conc_PO,conc_porucha, t);

//        if(t==(N_ITER-1)){
//            for(j=0;j<TYPY_PO;j++){
//                for(i=0;i<50;i++){
//                    conc_ve_smeru[0][j][i]=conc[j][150+i+150*Nx_UZLY].stara;
//                    conc_ve_smeru[1][j][i]=conc[j][150-i+150*Nx_UZLY].stara;
//                    conc_ve_smeru[2][j][i]=conc[j][150+(150+i)*Nx_UZLY].stara;
//                    conc_ve_smeru[3][j][i]=conc[j][150+(150-i)*Nx_UZLY].stara;
//                }
//            }
//            for(i=0;i<50;i++){
//                for(j=0;j<4;j++){
//                    fprintf(koncentrace_smer,"%e %e  ",conc_ve_smeru[j][0][i],conc_ve_smeru[j][1][i]);
//                }
//                fprintf(koncentrace_smer,"\n");
//            }
//            fclose(koncentrace_smer);
//        }
        /*
        for(n=0;n<N_PORUCHA;n++){
            x_porucha_average+=p_porucha[n].x[0];
        }
        x_porucha_average/=N_PORUCHA;
        */
        //average 1000cyklu rychlostni rozdeleni
//printf("things \n");
        for(k=0;k<TYPY_PO;k++){
            for(i=0;i<1000;i++){
                v_help_roz[k][i] = 0.0;
            }
        }
        for(n=0;n<N_PO;n++){
            if(p_oblast[n].v_cel<v_max[p_oblast[n].typ]){
                N_ofPartType[p_oblast[n].typ]++;
                v_help_roz[p_oblast[n].typ][(int)(p_oblast[n].v_cel/v_dilek[p_oblast[n].typ])]++;
            }
        }
        for(k=0;k<TYPY_PO;k++){
            for(i=0;i<1000;i++){
                v_help_roz[k][i]/=N_ofPartType[k];
            }
            for(i=0;i<1000;i++){
                v[k][i]=((t%IterBetweenSave)*v[k][i]+v_help_roz[k][i])/((t%IterBetweenSave)+1); // 100000 or 10000
            }
        }

        if(t%IterBetweenSave==0){ // 100000 or 10000
                //printf("jeje \n");
            /*iterace_help = t;
            testkdy = t/IterBetweenSave; // 100000 or 10000
            sprintf(testfile, "rozdeleni_PO_%d.txt", testkdy);
            if ((rychlost_test     =  fopen(testfile,"w"))        == NULL) { printf("Soubor rozdeleni_PO_%d.txt se nepodarilo otevrit. \n",testkdy);}
            for(i=0;i<1000;i++){
                //fprintf(rychlost_test, "%e %e %e %e\n", i*0.001*v_max[0], v[0][i], i*0.001*v_max[1], v[1][i]);
                for(k=0;k<TYPY_PO;k++){
                    fprintf(rychlost_test, "%i %e ", i, v[k][i]);
                }
                fprintf(rychlost_test, "\n");
            }
            fclose(rychlost_test);*/
            for(k=0;k<TYPY_PO;k++){
                for(i=0;i<1000;i++){
                    v[k][i] = v_help_roz[k][i];
                }
            }
        }
        for(k=0;k<TYPY_PO;k++){
            for(i=0;i<1000;i++){
                v_help_roz[k][i] = 0.0;
            }
        }
        if((t!=0)&&(t%IterBetweenSave==0)){ // 100000 or 10000
            //printf("jeje2 \n");
            for(k=0;k<TYPY_PO;k++){
                T_tot[k]=0.0;
            }
            for(n=0; n<N_PO; n++){
                T_tot[p_oblast[n].typ] += 0.5*(p_oblast[n].v_cel*p_oblast[n].v_cel*hmotnost[p_oblast[n].typ]);
            }
            for(k=0;k<TYPY_PO;k++){
                T_tot[k]/=(1.5*N_ofPartType[k]*kB);
            }
            //printf("T PO: %e, I: %e \n", T_tot[0], T_tot[1]);
            //fprintf(energie_PO,"%f %f \n",T_tot[0],T_tot[1]);

            PO_T_e=T_tot[0];
            PO_T_i=T_tot[1];
            fprintf(po_energie,"%e  %e  \n", PO_T_e,PO_T_i);

        }
        //fprintf(test_position,"%e %e %e %e %e %e %i\n", p_oblast[2].x[0],p_oblast[2].x[1],p_oblast[2].v[0],p_oblast[2].v[1],p_oblast[2].v[2],p_oblast[2].v_cel,p_oblast[2].typ);
        //fprintf(test_zdroj,"%e %e %e %e %e %e %i\n", p_zdroj_1[1].x[0], p_zdroj_1[1].x[1],p_zdroj_1[1].v[0],p_zdroj_1[1].v[1],p_zdroj_1[1].v[2],p_zdroj_1[1].v_cel,p_zdroj_1[1].typ);
//printf("pocet castic v oblasti %i ze zdroje %i\n", N_PO , N_out);
//if(t==5000)  {    save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, 4);  } /* 100 nano sekund  */
        //if(t==10000)  {    save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, 5);  } /* 250 nano sekund  */
        if(t%IterBetweenSave==(IterBetweenSave-1))  {   printf("save1 \n"); save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, (t+1)/IterBetweenSave);  } // 100000 or 10000 && 99999 or 9999
        if(t%IterBetweenSave==0)  {    printf("save2 \n"); save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, 2000+t/IterBetweenSave);  }
        //save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, 600+t);
//        if(t==20000) {    save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, 2);  } /* 1   mikro sekund */
//
//        if(t==30000){    save_n_E_U(p_U_old, p_E_old, conc,conc_PO,conc_porucha, t, 3);  } /* 2.5 mikro sekund */
        if((t%500000==0) && (t!=0))
        {
            printf("Do you want to stop computing y / n?\n");
            escape_char = getc(stdin);
            if( escape_char == escape_control) N_ITER=t;
        }
        t++;
    } while (t<N_ITER);

    /* konec casu */
 if(init_iterace<N_ITER) save_param(p_oblast,p_zdroj_1,p_zdroj_2,p_zdroj_3,p_zdroj_4,N_PO,N_ZDROJ,(N_ITER-1));
 /* vypocet proudu */
 for(i=0;i<POCET_SOND;i++){
    PROUD_CELKOVY[i]=0.0;
    for(j=0;j<TYPY_PO;j++){
        PROUD_CELKOVY[i]+=(PROUD_SONDY[j][i]*naboj[j])/((N_ITER-iterace_ustaleni)*DT[j]*(2 * PI * R_S[i] * Z_PO));
    }

 }
 /*
 for(k=0;k<TYPY_PO;k++){
    for(i=0;i<POCET_SOND;i++){
        for(j=0;j<DELENI_SONDA;j++){
            TOK_SONDY[k][i][j]=((POMOC_TOK[k][i][j]*naboj[k])/((N_ITER-iterace_ustaleni)*DT[k]*((2*PI*R_S[i]*Z_PO)/DELENI_SONDA)));
        }
    }
 }
 for(j=0;j<DELENI_SONDA;j++){
    for(i=0;i<POCET_SOND;i++){
        for(k=0;k<TYPY_PO;k++){
            fprintf(tok_sonda,"%e ",TOK_SONDY[k][i][j]);
        }
    }
    fprintf(tok_sonda,"\n");
 }*/
 fclose(tok_sonda);

 clock_t end=clock();
 double time_in_seconds = (double) (end - start)/CLOCKS_PER_SEC;

 for(k=0;k<TYPY_PO;k++){
    N_ofPartType[k] = 0;
 }
 for (n=0;n<N_PO;n++){
    N_ofPartType[p_oblast[n].typ]++;
 }

 for(k=0;k<TYPY_PO;k++){
    printf("particles of part type %i: %i \n", k, N_ofPartType[k]);
 }


 cas_p=cas_p/N_ITER; cas_k=cas_k/N_ITER; cas_r=cas_r/N_ITER; cas_s=cas_s/N_ITER; cas_z=cas_z/N_ITER;

 dny= (int) (floor(time_in_seconds/86400));
 hodiny= (int) (floor((fmod(time_in_seconds,86400))/3600));
 minuty= (int) (floor(fmod(fmod(time_in_seconds,86400),3600)/60));
 sekundy= (int) fmod(fmod(fmod(time_in_seconds,86400),3600),60);

 fprintf(soubor,"program trval dni: %i, hodin: %i, minut: %i, sekund: %i \n", dny, hodiny, minuty, sekundy);
 fprintf(soubor,"v pracovni oblasti bylo na konec: %i castic \n", N_PO);
  for(k=0;k<TYPY_PO;k++){
    fprintf(soubor,"No of %i type: %i ",k, N_ofPartType[k]);
 }
 fprintf(soubor,"\n");
 fprintf(soubor,"LambdaD %e ",Lambda_D);
 Lambda_D = sqrt((eps0*kB/pow(naboj[0],2))/((n_NENARUS[0]/T_tot[0])+(n_NENARUS[1]/T_tot[1])));
 fprintf(soubor,"with final T LambdaD %e \n",Lambda_D);

 for(i=0;i<POCET_SOND;i++){
    fprintf(soubor,"napeti_%i %f \n",i,U_S[i]);
    for(j=0;j<TYPY_PO;j++){
        fprintf(soubor,"proud_%i_typu %i \n",j,PROUD_SONDY[j][i]);
    }
    fprintf(soubor,"celk_proud %.9e \n",PROUD_CELKOVY[i]);
    fprintf(soubor,"\n");
 }
 fprintf(soubor,"polohy: %f \nkoncentrace: %f \nrychlosti: %f \nsila: %f \nzdroj: %f \n", cas_p, cas_k, cas_r, cas_s, cas_z);
 fclose(soubor);
// fclose(energie_PO);
fclose(kontrola);
fclose(vyvoj_proudu);
fclose(zdroj_energy);
fclose(zdroj_tok);
fclose(po_tok);
fclose(po_energie);

fclose(test_position);
fclose(test_zdroj);
fclose(vyvoj_naboje);
umfpack_di_free_numeric(&Numeric);

printf("Uvolnovani pameti ...\n");
for (i=0; i<TYPY_PO; i++) {free((void **) conc[i]);}
free((void *) conc);
for (i=0; i<counter_e_ela; i++) {free((void **) sigma_e_ela[i]);}
free((void *) sigma_e_ela);
for (i=0; i<counter_e_exc; i++) {free((void **) sigma_e_exc[i]);}
free((void *) sigma_e_exc);
for (i=0; i<counter_e_ion; i++) {free((void **) sigma_e_ion[i]);}
free((void *) sigma_e_ion);
free((void *) p_zdroj_1);
free((void *) p_zdroj_2);
free((void *) p_zdroj_3);
free((void *) p_zdroj_4);
free((void *) p_oblast);
free((void *) p_porucha);
free((void *) p_flag);
free((void *) p_constpot);
free((void *) p_constpot_init);
free((void *) p_U);
free((void *) p_gradU_BC);
free((void *) p_U_Bound_pointQ);
free((void *) p_U_PIC);
free((void *) p_U_old);
free((void *) p_E);
free((void *) p_E_old);
free((void *) vypadle_ind);
free((void *) p_Ax);
free((void *) p_Ai);
free((void *) p_Ap);
free((void *) p_b);

printf("Konec programu 2D modelu. \n");
return 0;
}
