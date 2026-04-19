#include <stdlib.h>
#include <stdio.h>
#include <pcg_variants.h>
#include "global.h"


/* definice globalnich promennych */
/* elektrony typ castice 0, ionty typ castice 1 */ /* 10% elneg*/
int    POCET_PO[TYPY_PO] 		= {2000000,		2000000};  /* pocet castic v PO 	puvodni je 2M pro kazdy	.. pro vetsi je 8M			*/
int    POCET_ZDROJ[TYPY_ZDROJ]  = {15000,		15000  };	/* pocet castic ve ZDROJI, mam 4 pole s tímto poctem castic - zdroj pres x a zdroj pres y -- puvodni je 15k	pro vetsi 60K*/
double n_NENARUS[TYPY_PO]		= {1.59e15, 	1.59e15}; 	/* koncentrace daneho druhu v nenarusenem plazmatu - kvuli inicializaci PO 	9.54e14, 	1.59e15, 6.36e14		*/

int    N_ITER 			        = 1000001;                                                         /* pocet casovych iteraci */ /* dynamic 2000001 na 311 */
int    iterace_ustaleni         = 10000000;

double TEPLOTA_PO[TYPY_PO]		= {33500.0,			300.0};	/* teplota v PO ze stredni energie 4.27eV = 33034.27K  (pro 100V/m pri tlaku 133Pa =3.11Td .. <E>=3/2 KbT)					*/
double TEPLOTA_ZDROJE[TYPY_PO]	= {33500.0,			300.0};	/* teplota ve ZDROJI 					*/

double hmotnost[TYPY_PO]		= {9.10938215e-31, 	6.6353628e-26};	/* hmotnost castic Ar+ 6.6353628e-26, Ne+_like	3.3206524e-26, Kr+_like 1.39149851e-25 */
double DT[TYPY_PO]				= {5e-12, 		5e-12};	/* casovy krok         1e-11 ion on 311 else 1e-8        */
double DT_zdroj[TYPY_PO]		= {5e-12, 	    5e-12};	/* casovy krok ve zdroji        */
double naboj[TYPY_PO]			= {-1.60217657e-19, 1.60217657e-19};	/* naboj castic		            */

/* hrebinek sond, sond je 2 x 3 */

const int SONDA[POCET_SOND]     = { 1};//,       1,       1};//,       1,       1,       1,       1,       1,       1      };
const double R_S[POCET_SOND]    = { 1e-4};//, 3.01e-4, 3.01e-4};//, 3.01e-4, 3.01e-4, 3.01e-4, 3.01e-4, 3.01e-4, 3.01e-4};//1e-4
const double X_S[POCET_SOND]    = { 0.0};//,    0.0,     0.01};//,   0.0,     0.0,     0.0,     0.008,   0.008,   0.008    };
const double Y_S[POCET_SOND]    = { 0.0};//,     0.0,     0.0};//,  -0.008,   0.0,     0.008,  -0.008,   0.0,     0.008   };
      double U_S[POCET_SOND]    = { 6.0};//,     0.0,     0.0};//,   0.0,     0.0,     0.0,    -0.008,  -0.008,  -0.008    }; // 6V na 311


/* double MAX_RYCH[TYPY_PO]		= {3.0e6, 			3.0e4			};	 max.rychlost v histogramech 	*/
/* int	   SRAZKY[TYPY_PO]		= {3,				2				};   pocet srazkovych procesu pro dany typ castic, pro vypis do souboru frekvence 	*/
/* int	   MAX_SRAZKY			=  3;									 max pocet typu srazek pro ! druh castic 										*/
//double lambda_MAX[TYPY_PO]    = {1.95E-4,         3.279226E-5};//,    3.279226E-5};	/* celkova stredni volna draha nabitych castic, 0 .. elektron, 1 .. iont    		*/
double M_MAX[TYPY_PO]		    = {3.5e-13, 		3.8e-15}; 	/* jednotka m3/s, max{rychlost*sigma}	*/
                                                    /* M_MAX_el  = v*sigma for 12,1eV ..... IT CAN BE MODIFIED IF MAIN!! based on actual crosssections used*/
													/* M_MAX_IO  = V_MAX_IO * (P_PRENOS_IO + P_ROZPTYL_IO) 		*/
													/* M_MAX_IO  = 4000.0 * (5.5e-19 + 4.0e-19) = 3.8e-15 		*/
													/* M_MAX_NEG = V_MAX_NEG * (P_ELASTIC_NEG + nulova srazka) 	*/
													/* M_MAX_NEG = 4000.0 * (4.3e-19 + 5.2e-19) = 3.8e-15 		*/
/* double P_PRENOS_IO 			= 5.5e-19;		*/	/* uc.prurez charge transferu, pro klad. ionty 				*/
/* double P_ELASTIC_NEG			= 4.3e-19;		*/	/* uc.prurez elastickeho rozptylu, pro neg.ionty 			*/

/* outputy */
FILE           *potencial_w, *intenzita_w, *koncentrace_w, *koncentrace_PO, *koncentrace_porucha, *rychlost_PO, *velocity_source, *rychlost_porucha, *radial_velocity, *angular_velocity, *zdroj_energy,*rychlostx_PO,*rychlosty_PO,*rychlostz_PO;
FILE           *kontrola, *oblast, *zdroj_1, *zdroj_2, *zdroj_3, *zdroj_4, *parametry;
FILE           *oblast_1, *zdroj_1_1, *zdroj_2_1, *zdroj_3_1, *zdroj_4_1, *parametry_1;
FILE           *koncentrace_smer, *tok_sonda, *zdroj_tok, *po_tok, *po_energie;

DDCASTICE      *p_zdroj_1, *p_zdroj_2, *p_zdroj_3, *p_zdroj_4, *p_oblast, *p_porucha, *pAllSourcesSteadyVelocity;
double 	       *p_U,        *p_U_old, *p_gradU_BC, *p_U_Bound_pointQ;
double		   *p_U_PIC;
double         **sigma_e_ela, **sigma_e_exc, **sigma_e_ion;
INTENZITA      *p_E,        *p_E_old;
KONCENTRACE    **conc, **conc_PO, **conc_porucha;   	/* dvourozmerne pole pro ukladani objemove koncentrace */
double         conc_ve_smeru[4][TYPY_PO][50];
//double         **conc_DD;								/* 2D pole pro ukladani radialni koncentrace */
//double 		   **energie_PO;
//VEKTOR    	   **rychlosti_PO;
//VEKTOR 		   **toky_PO;
//double		   ***frekvence_PO, ***frekvence_PO_OLD;	/* srazkove frekvence na siti */
//double		   **hist_PO, **hist_SHEATH, **hist_OUT;
INTENZITA      E_extern;
VEKTOR         a_zdroj[TYPY_PO];
int    t; /* cas - iterace */
int    average_counter;
int    N_PO, N_ZDROJ, N_PORUCHA;
int    POCET_PORUCHA[TYPY_PO];
int	   N_aktual[TYPY_PO];
double    Scattered[TYPY_PO];
double    Scattered_zero;
double    Lambda_D;
double    Omega_p[TYPY_PO];
double V_kvadr[TYPY_PO];
double PartForCellInBulk[TYPY_PO];
double Qmdiv[TYPY_PO];
double KONSTANTA_POHYBU[TYPY_PO];
double vmax[TYPY_PO];
double v_radial[TYPY_PO][40][30],v_angular[TYPY_PO][40][30], v_distribution[TYPY_PO][1000],vx_distribution[TYPY_PO][1000],vy_distribution[TYPY_PO][1000],vz_distribution[TYPY_PO][1000];

double Omega_c[TYPY_PO];
double Brot_11[TYPY_PO],Brot_12[TYPY_PO],Brot_13[TYPY_PO],Brot_21[TYPY_PO],Brot_22[TYPY_PO],Brot_23[TYPY_PO],Brot_31[TYPY_PO],Brot_32[TYPY_PO],Brot_33[TYPY_PO];

double TAU_MAX[TYPY_PO];

/* mala_sonda_1	0, mala_sonda 1, rovinna sonda 2 */
int 	POMOC_TOK[TYPY_PO][POCET_SOND][DELENI_SONDA];
double  TOK_SONDY[TYPY_PO][POCET_SOND][DELENI_SONDA];
//double  VRMS_SONDY[TYPY_PO][3][DELENI_MAX];
//double	ENERGIE_SONDY[TYPY_PO][3][DELENI_MAX];
//double  HISTOGRAMY_SONDY[TYPY_PO][3][POCET_SLOUPCU_SONDY];
int 	PROUD_SONDY[TYPY_PO][POCET_SOND];
int     POMOC_PROUD[TYPY_PO][POCET_SOND][20];
double  ARITM_PROUD[TYPY_PO][POCET_SOND];
double  ARITM_ODCH_PROUD[TYPY_PO][POCET_SOND];
double 	PROUD_CELKOVY[POCET_SOND];

/* ostatni promenne */

double Z_PO;
double h_PO;
double Cell_volume;

int    Ny_UZLY;
int    Ny_UZLY_PIC;

double X0_PO;
double Y0_PO;

double Z_i, Z_f;
double Y_ZDROJ_1, Y_ZDROJ_2, Y_ZDROJ_3, Y_ZDROJ_4;
double X_ZDROJ_1, X_ZDROJ_2, X_ZDROJ_3, X_ZDROJ_4;

double Y_PORUCHA;
double Z_PORUCHA;

double Tot_charge_in_space, Real_Qpoint_value;

double sigma_e_tot;

int counter_e_ela, counter_e_exc, counter_e_ion;

int    N_out;
int    N_ROW;
int    N_NOZERO;
int    N_ROW_PIC;
int    N_NOZERO_PIC;
int    N_REZERVA_OBLAST;
int    N_REZERVA_VYPADLE;
int    N_BC_NODES;

double Rand_modif;
pcg64_random_t rng[NT];
int Nr_thread;




