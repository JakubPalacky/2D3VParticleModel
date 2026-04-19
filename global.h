#include <stdlib.h>
#include <stdio.h>
#include <pcg_variants.h>

#ifndef GLOBAL_H                                                                /* podmineny preklad proti opakovanemu inkludovani */
#define GLOBAL_H
#endif

/* vstupni parametry */
#define n0 					1.59e15     /* 1.59e15 pri 133 Pa, 8.73e14 pri 73 Pa, 0.91e15 pri 67 Pa, 2.00e15 pri 270 Pa , */ /* koncentrace nenaruseneho plazmatu */
#define n_NEUT 				3.21e22		/* 3.21e22 pri 133 Pa, 1.61e22 pri 66.5Pa 4.82e22 pri 199.5Pa..9.63e22 pri 399,97Pa, 1.61e23 pøi 665Pa, 3.21e20 pøi 1,33Pa              1.76e22 pri 73 Pa, 1.62e22 pri 67 Pa, 6.52e22 pri 270 Pa*/ /* koncentrace neutralu (rozptylovych center) */
#define Nx_UZLY  			401                                                 /* pocet uzlu na x-ove hrane PO .. pocet uzlu +2 */
#define X_PO 		 		0.01                                                /* hrana x obdelnikove pracovni oblasti */
#define Y_PO 				0.01                                                /* hrana y obdelnikove pracovni oblasti */
                                                                               /* !!!! Y_PO MUSI BYT NASOBEK h_PO !!!! */
#define TYPY_PO 			2                                                   /* pocet typu castic v PO */
#define TYPY_ZDROJ 			2                                                   /* typ 0 ... elektrony, typ 1 ... ionty */
#define POCET_SOND          1                                                   /* pocet valcovych sond v PO */

#define EL_POLE             100.0                                                 /* vnejsi el. pole [V/m] bude pevne ve smeru x ... pro 66,5Pa tedy 0,5 torr pouzit 50 a pro 1,5torr tedy 199,5Pa pouzit 150*/

/* MAG. Field in model */
#define MAGNETIC_FIELD      0                                                   /* mag. field off 0 on 1 */
#define BFIELD_Z            0                                               /* mag. field in z direction [T]*/
#define BFIELD_X            0.01                                                /*mag. field in x direction... for mag sheath*/

#ifdef _OPENMP
  #define NT				30												    /* pocet procesoru */
  #else
    #define NT              1
#endif

#define IterBetweenSave     10000
																				/* param. male_sondy_1 */
#define PERIOD				0		                                            /* zda jsou v modelu periodicke hranice */

/* parametry rovinne poruchy */
#define PORUCHA             0           /* 1 .. porucha on, 0 .. off */
#define T_PORUCHA           20000         /* casovy krok kdy vstoupi porucha do PO */
#define n_porucha           1.59e15     /* porucha bude jednou tolik */
#define X_PORUCHA           0.005       /* sirka poruchy na 10 bunek */


/* inicializace nebo cteni ze souboru */

#define CTENI_ZE_SOUBORU 	1   	  		    								    /* zda se ma p_oblast nacist ze souboru */

/* postprocessing */

#define KROK_PRUMEROVANI 	10 		 		    /* vratit 10 */			    		    /* s jakym krokem se prumeruje */
#define PRUMEROVANI 		0       	    /* vratit 100000 */		    	     	/* od ktere iterace zacina prumerovani vysledku */
//#define RESET_PRUMEROVANI   1000
/* radialni koncentrace */

#define N_R 				800  											    	/* pocet intervalu vrstvy pro postprocessing koncentrace */

/* vypocet proudu a toku na sondy */

/* #define VYPOCET_PROUDU 		500 		*/										/* interval pro vypocet proudu a toku */
/* #define TOK_START			170000		*/										/* kterou iteraci zacina vypocet uhlove zavislosti toku */
#define DELENI_SONDA			360													/* pocet dilku na obvodu sondy */
/* #define DELENI_MAX			400			*/										/* max pocet dilku na povrsich sond */
/* #define POCET_SLOUPCU_SONDY 80 */

/* histogramy */
#define Initial_source_time 1e-7 //change back to -7

/* #define POCET_SLOUPCU       80	*/												/* pocet sloupcu histogramu */
/* #define R_SHEATH			0.003	*/											/* hranice sheathu a PO */
/* #define R_OUT				0.02	*/											/* vnejsi hranice oblasti, kde se pocita historgram */

/* nastaveni vypoctu */
/* #define KROK_UKLADANI 		1000   */          	    							/* s jakym krokem se ukladaji pole koncentrace a potencial */
/* #define FREKVENCE_ITERACE	100		*/											/* pocet iteraci pred PRUMEROVANI, kdy se zacnou pocitat srazkove frekvence na siti */
#define N_PRIDAVEK_OBLAST 	100000                                               /* pocet pridanych castic pri preteceni p_oblast */
#define N_PRIDAVEK_VYPADLE 	5000    											/* pocet pridanych castic pri preteceni vypadle_ind */

/* fyzikalni konstanty */
#define PI 					3.14159265358979
#define kB 					1.380658e-23
#define eps0 				8.854187817e-12
#define M_Ar_au             40
#define M_Used_au           40 /* Neon like 20 , Kr like 84*/

/* prevodi konstanty */
#define JtoEV               6.2415091e18
#define EVtoJ               1.602176634e-19

/* definice globalnich typu */
typedef struct{                                                                 /* definice castice */
        double x[3];
        double v[3];
        double v_cel;
        double t_s;                                                             // time to scatter
        int typ;                                                                /* typ 0 ... elektrony, typ 1 ... ionty */
        int in;                                                                 /* vstoupila uz castice poruchy do PO 0..false 1..true */
        } DDCASTICE;

typedef struct{																	/* definice pole intenzity E */
        double x;
        double y;
        } INTENZITA;

typedef struct{
		double x;
		double y;
		double z;
		} VEKTOR;

typedef struct{
        double new;
        double stara;
        } KONCENTRACE;


/* deklarace globalnich promennych */
extern  int    POCET_PO[TYPY_PO];
extern  int    POCET_ZDROJ[TYPY_PO];
extern  double TEPLOTA_PO[TYPY_PO];
extern  double TEPLOTA_ZDROJE[TYPY_PO];
extern  double hmotnost[TYPY_PO];
extern  double DT[TYPY_PO];
extern  double DT_zdroj[TYPY_PO];
extern  double naboj[TYPY_PO];
extern  double vmax[TYPY_PO];
extern  double v_radial[TYPY_PO][40][30],v_angular[TYPY_PO][40][30], v_distribution[TYPY_PO][1000],vx_distribution[TYPY_PO][1000],vy_distribution[TYPY_PO][1000],vz_distribution[TYPY_PO][1000];
extern  double Omega_c[TYPY_PO];
extern  double Brot_11[TYPY_PO],Brot_12[TYPY_PO],Brot_13[TYPY_PO],Brot_21[TYPY_PO],Brot_22[TYPY_PO],Brot_23[TYPY_PO],Brot_31[TYPY_PO],Brot_32[TYPY_PO],Brot_33[TYPY_PO];
extern  double TAU_MAX[TYPY_PO];
extern double    Scattered[TYPY_PO];
extern double    Scattered_zero;

extern const int SONDA[POCET_SOND];
extern const double R_S[POCET_SOND];
extern const double X_S[POCET_SOND];
extern const double Y_S[POCET_SOND];
extern double U_S[POCET_SOND];

extern  double n_NENARUS[TYPY_PO];
//extern  int	 SRAZKY[];
//extern  int	 MAX_SRAZKY;
//extern  double MAX_RYCH[];
//extern  double lambda_MAX[TYPY_PO];
extern  double M_MAX[TYPY_PO];
//extern  double P_PRENOS_IO;
//extern  double P_ELASTIC_NEG;

extern FILE            *potencial_w, *intenzita_w, *koncentrace_w, *koncentrace_PO, *koncentrace_porucha, *rychlost_PO, *velocity_source, *rychlost_porucha, *radial_velocity, *angular_velocity, *zdroj_energy, *rychlostx_PO,*rychlosty_PO,*rychlostz_PO;
FILE           *kontrola, *oblast, *zdroj_1, *zdroj_2, *zdroj_3, *zdroj_4, *parametry;
extern FILE            *kontrola, *oblast, *zdroj_1, *zdroj_2, *zdroj_3, *zdroj_4, *parametry;
extern FILE            *oblast_1, *zdroj_1_1, *zdroj_2_1, *zdroj_3_1, *zdroj_4_1, *parametry_1;
extern FILE            *koncentrace_smer, *tok_sonda, *zdroj_tok, *po_tok,*po_energie;

extern  DDCASTICE      *p_zdroj_1, *p_zdroj_2, *p_zdroj_3, *p_zdroj_4, *p_oblast, *p_porucha, *pAllSourcesSteadyVelocity;
extern  double 	       *p_U,        *p_U_old, *p_gradU_BC,*p_U_Bound_pointQ;
extern  double		   *p_U_PIC;
extern  double         **sigma_e_ela, **sigma_e_exc, **sigma_e_ion;
extern  INTENZITA      *p_E,        *p_E_old;
extern  KONCENTRACE    **conc, **conc_PO, **conc_porucha;
extern  double         conc_ve_smeru[4][TYPY_PO][50];
//extern  double         **conc_DD;
//extern  double 		   **energie_PO;
//extern  VEKTOR    	   **rychlosti_PO;
//extern  VEKTOR 		   **toky_PO;
//extern  double		   ***frekvence_PO, ***frekvence_PO_OLD;
//extern  double		   **hist_PO, **hist_SHEATH, **hist_OUT;
extern  INTENZITA      E_extern;
extern  VEKTOR         a_zdroj[TYPY_PO]; // zrchleni ve zdroji s extern E field
extern  int            t, N_ITER,iterace_ustaleni;
extern  int            average_counter;
extern  int 		   N_PO, N_ZDROJ, N_PORUCHA;
extern  int            POCET_PORUCHA[TYPY_PO];
extern  int			   N_aktual[TYPY_PO];
extern  double		   V_kvadr[TYPY_PO];
extern  double         PartForCellInBulk[TYPY_PO];
extern  double         Qmdiv[TYPY_PO];
extern  double         KONSTANTA_POHYBU[TYPY_PO];
extern double    Lambda_D;
extern double    Omega_p[TYPY_PO];

extern  int 	    POMOC_TOK[TYPY_PO][POCET_SOND][DELENI_SONDA];
extern  double	    TOK_SONDY[TYPY_PO][POCET_SOND][DELENI_SONDA];
//extern  double    VRMS_SONDY[TYPY_PO][3][DELENI_MAX];
//extern  double	ENERGIE_SONDY[TYPY_PO][3][DELENI_MAX];
//extern  double    HISTOGRAMY_SONDY[TYPY_PO][3][POCET_SLOUPCU_SONDY];
extern  int 	PROUD_SONDY[TYPY_PO][POCET_SOND];
extern  int     POMOC_PROUD[TYPY_PO][POCET_SOND][20];
extern  double  ARITM_PROUD[TYPY_PO][POCET_SOND];
extern  double  ARITM_ODCH_PROUD[TYPY_PO][POCET_SOND];
extern  double 	PROUD_CELKOVY[POCET_SOND];

extern  double  Z_PO;
extern  double  h_PO;
extern  double  Cell_volume;

extern  int     Ny_UZLY;
extern  int     Ny_UZLY_PIC;

extern  double  X0_PO;
extern  double  Y0_PO;

extern  double Z_i, Z_f;
extern  double Y_ZDROJ_1, Y_ZDROJ_2, Y_ZDROJ_3, Y_ZDROJ_4;
extern  double X_ZDROJ_1, X_ZDROJ_2, X_ZDROJ_3, X_ZDROJ_4;

extern  double  Y_PORUCHA;
extern  double  Z_PORUCHA;

extern double Tot_charge_in_space, Real_Qpoint_value;
extern int counter_e_ela, counter_e_exc, counter_e_ion;

extern  int     N_out;
extern  int     N_ROW;
extern  int     N_NOZERO;
extern  int     N_ROW_PIC;
extern  int     N_NOZERO_PIC;
extern  int     N_REZERVA_OBLAST;
extern  int     N_REZERVA_VYPADLE;
extern  int     N_BC_NODES;

extern  double Rand_modif;
extern  pcg64_random_t rng[NT];
extern  int Nr_thread;


