#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "PIC.h"

void FLAG_CONSTPOT(int *prip_flag, double *prip_constpot)
{ int i, j, k, l;
  int i_0, j_0, i_1, j_1;
  double x_uzel, y_uzel, distance;
  double radius = X_PO/2;

    /* flag = 0 - uzel ma konst.potencial, 1 - ma 4 sousedy,                                */
    /*  	    2+i (i=0, i<POCET_SOND) - okoli i-te sondz             						  */
    for (j=0; j<Ny_UZLY; j++){
        for(i=0;i<Nx_UZLY;i++){
            k=i+j*Nx_UZLY;
            prip_flag[k] 	 = 1;                          	/* uzel ma vsechny ctyri sousedy */
            prip_constpot[k] = 0.0;
        }
    }
    /* okrajove podminky - el pole */
    for (k = 0; k < Nx_UZLY; k++) {
        prip_flag[k] 						 			 = 0;
        prip_flag[(Ny_UZLY_PIC-1)*Nx_UZLY+k] 			 = 0;
        prip_constpot[k]                                 +=EL_POLE*((k*h_PO)+X0_PO);
        prip_constpot[(Ny_UZLY_PIC-1)*Nx_UZLY+k]         +=EL_POLE*((k*h_PO)+X0_PO);
        //printf("%f \n",prip_constpot[k]);
        //printf("%i %f \n",((Ny_UZLY_PIC-1)*Nx_UZLY+k), prip_constpot[(Ny_UZLY_PIC-1)*Nx_UZLY+k]);
    }
    for (k = 1; k < Ny_UZLY_PIC-1; k++) {
        prip_flag[k*Nx_UZLY] 	                         = 0;
        prip_flag[(k+1)*Nx_UZLY-1]                       = 0;
        prip_constpot[k*Nx_UZLY] 	                     += EL_POLE*(X0_PO);
        prip_constpot[(k+1)*Nx_UZLY-1]                   += EL_POLE*(X0_PO+((Nx_UZLY-1)*h_PO));
        //printf("%f \n",prip_constpot[k*Nx_UZLY]);
        //printf("%f \n",prip_constpot[(k+1)*Nx_UZLY-1]);
    }

    /* okoli vice valcovych sond */
    for(l=0; l<POCET_SOND; l++){
        if (SONDA[l]==1) {
            i_0 = floor((X_S[l]-R_S[l]-X0_PO )/h_PO);    i_1 = ceil((X_S[l]+R_S[l]-X0_PO )/h_PO);
            j_0 = floor((Y_S[l]-R_S[l]-Y0_PO )/h_PO);    j_1 = ceil((Y_S[l]+R_S[l]-Y0_PO )/h_PO);

            for (i = i_0; i <= i_1; i++) {
                for (j = j_0; j <= j_1; j++) {
                    k = i + Nx_UZLY*j;
                    x_uzel = X0_PO + i * h_PO;
                    y_uzel = Y0_PO + j * h_PO;
                    if (sqrt(pow(X_S[l]-x_uzel,2)+pow(Y_S[l]-y_uzel,2)) <= R_S[l]) {
                        prip_flag[k] = 0;
                        prip_constpot[k] += U_S[l];
                    }
                    else {  prip_flag[k] = l+2;  }
                }
            }
        }
    }
}

void CONSTPOT_BOUNDARY_CHANGE(double *prip_constpot_init, double *pole_constpot, double *pole_U_BC)
{
    int k;
    for(k=0;k<N_ROW;k++){
        //prip_constpot_init[k] = pole_constpot[k];
        pole_constpot[k] = prip_constpot_init[k] - pole_U_BC[k];
    }
}


void PRIPRAV_B(DDCASTICE *prip_PO, DDCASTICE *prip_porucha, double *prip_b, double *prip_constpot, KONCENTRACE **prip_conc,KONCENTRACE **prip_conc_PO, KONCENTRACE **prip_conc_porucha)      /* pole priprav_PO obsahuje castice v PO */
{                                                                               						/* pole priprav_b ma (Nx_UZLY * Ny_UZLY) polozek typu double */
  int  		i,j,k,l,m,n;
  double 	i_notrunc,j_notrunc;
  double 	dx_h,dy_h;
  double 	celkovy_naboj, real_cell_volume;//, objem_bunky;

  for (l = 0; l < TYPY_PO; l++) {												/* vynulovani koncentrace */
    for (k = 0; k < N_ROW; k++) {
        prip_conc[l][k].new     = 0.0;
        prip_conc_PO[l][k].new       = 0.0;
        prip_conc_porucha[l][k].new  = 0.0;
    }
  }

    KONCENTRACE **prip_conc0;

    if ( (prip_conc0 = (KONCENTRACE **) 		malloc(TYPY_PO*sizeof(KONCENTRACE *)) ) == NULL ) 	{ printf("Malo pameti pro pole conc."); 			exit(1); };

   	for (m = 0; m < TYPY_PO; m++) {
   	  if ( (prip_conc0[m] = (KONCENTRACE *) 	malloc(N_ROW*sizeof(KONCENTRACE)) ) 	== NULL ) 	{ printf("Malo pameti pro %d.radek pole conc.", m); exit(1); };
   	}

   	for (m = 0; m < TYPY_PO; m++) {
   	  for (n = 0; n < N_ROW; n++) {
   	    prip_conc0[m][n].new = 0.0;
   	  }
   	}
    // border concentration fix: k = i+j*Nx_UZLY
//    for (m = 0; m < TYPY_PO; m++) {
//    	// right on edge sides must have +1/2*PartForCellInBulk and corner +3/4*PartForCellInBulk
//            //bottom j = 1 and top j = Ny_UZLY-2
//    	for (i = 1; i < Nx_UZLY-1; i++){
//            k = i;
//            prip_conc0[m][k].new += 0.5*PartForCellInBulk[m];
//            k = i+(Ny_UZLY-1)*Nx_UZLY;
//            prip_conc0[m][k].new += 0.5*PartForCellInBulk[m];
//    	}
//            //sides left i = 1 and right i = Nx_UZLY-2
//    	for(j = 1; j < Ny_UZLY-1; j++){
//            k = j*Nx_UZLY;
//            prip_conc0[m][k].new += 0.5*PartForCellInBulk[m];
//            k = (j+1)*Nx_UZLY-1;
//            prip_conc0[m][k].new += 0.5*PartForCellInBulk[m];
//    	}
//        //corners
//        //bottom left 0 , 0
//        k = 0;
//        prip_conc0[m][k].new += 0.75*PartForCellInBulk[m];
//        //bottom right Nx_UZLY-1 , 0
//        k = Nx_UZLY-1;
//        prip_conc0[m][k].new += 0.75*PartForCellInBulk[m];
//        //top left 0 , Ny_UZLY-1
//        k = (Ny_UZLY-1)*Nx_UZLY;
//        prip_conc0[m][k].new += 0.75*PartForCellInBulk[m];
//        //top right Nx_UZLY-1 , Ny_UZLY-1
//        k = (Ny_UZLY*Nx_UZLY)-1;
//        prip_conc0[m][k].new += 0.75*PartForCellInBulk[m];
//    }

   	#pragma omp parallel default(shared) private(k,i,j,dx_h,dy_h,i_notrunc,j_notrunc) num_threads(NT)
    {
        #pragma omp for
        for (l = 0; l < N_PO; l++) {                                                   /* CIC */
        k = 0;
          i_notrunc = (prip_PO[l].x[0]-X0_PO)/h_PO;                                    /* nascitavani castic pro jednotlive uzly */
          j_notrunc = (prip_PO[l].x[1]-Y0_PO)/h_PO;
          i = (int) floor(i_notrunc);
          j = (int) floor(j_notrunc);
          k = i + Nx_UZLY * j;
          // if ( (k + 1 + Nx_UZLY) >= N_ROW) printf("\n  Index (k + 1 + Nx_UZLY) ve fci PRIPRAV_B mimo mez N_ROW, k + 1 + Nx_UZLY = %d, N_ROW = %d \n", k + 1 + Nx_UZLY, N_ROW);
          dx_h = i_notrunc - i;
          dy_h = j_notrunc - j;
            if ((i<(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
                prip_conc0[prip_PO[l].typ][k+1].new         	+=  (dx_h)   * (1-dy_h);
                prip_conc0[prip_PO[l].typ][k+Nx_UZLY].new   	+=  (1-dx_h) * (dy_h);
                prip_conc0[prip_PO[l].typ][k+1+Nx_UZLY].new	    +=  (dx_h)   * (dy_h);
            }
            else if ((i==(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
                prip_conc0[prip_PO[l].typ][k+Nx_UZLY].new   	+=  (1-dx_h) * (dy_h);
            }
            else if ((i<(Nx_UZLY-1)) && (j==(Ny_UZLY-1))) {
                prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
                prip_conc0[prip_PO[l].typ][k+1].new         	+=  (dx_h)   * (1-dy_h);
            }
            else {
                prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
            }
        }
    }

    for (l = 0; l < N_ROW_PIC; l++) {
        for (k = 0; k < TYPY_PO; k++) {
            prip_conc_PO[k][l].new=prip_conc0[k][l].new;
        }
    }

    /* nacitani castic z poruchy do koncentrace */
    /*
    if((PORUCHA)&&(t>=T_PORUCHA)){
        #pragma omp parallel default(shared) private(k,i,j,dx_h,dy_h,i_notrunc,j_notrunc) num_threads(NT)
        {
            #pragma omp for
            for (l = 0; l < N_PORUCHA; l++) {
                if(prip_porucha[l].in==1){
                    k = 0;
                    i_notrunc = (prip_porucha[l].x[0]-X0_PO + h_PO)/h_PO;
                    j_notrunc = (prip_porucha[l].x[1]-Y0_PO + h_PO)/h_PO;
                    i = (int) floor(i_notrunc);
                    j = (int) floor(j_notrunc);
                    k = i + Nx_UZLY * j;
                    // if ( (k + 1 + Nx_UZLY) >= N_ROW) printf("\n  Index (k + 1 + Nx_UZLY) ve fci PRIPRAV_B mimo mez N_ROW, k + 1 + Nx_UZLY = %d, N_ROW = %d \n", k + 1 + Nx_UZLY, N_ROW);
                    dx_h = i_notrunc - i;
                    dy_h = j_notrunc - j;
//                    if ((i<(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
                        prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
                        prip_conc0[prip_PO[l].typ][k+1].new         	+=  (dx_h)   * (1-dy_h);
                        prip_conc0[prip_PO[l].typ][k+Nx_UZLY].new   	+=  (1-dx_h) * (dy_h);
                        prip_conc0[prip_PO[l].typ][k+1+Nx_UZLY].new	    +=  (dx_h)   * (dy_h);
//                    }
//                    else if ((i==(Nx_UZLY-1)) && (j<(Ny_UZLY-1))) {
//                        prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
//                        prip_conc0[prip_PO[l].typ][k+Nx_UZLY].new   	+=  (1-dx_h) * (dy_h);
//                    }
//                    else if ((i<(Nx_UZLY-1)) && (j==(Ny_UZLY-1))) {
//                        prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
//                        prip_conc0[prip_PO[l].typ][k+1].new         	+=  (dx_h)   * (1-dy_h);
//                    }
//                    else {
//                        prip_conc0[prip_PO[l].typ][k].new				+=  (1-dx_h) * (1-dy_h);
//                    }
                }
                else continue;
            }
        }
    }
    */
    /* konci paralelni cyklus nacitani koncentrace*/

  	  for (m = 0; m < TYPY_PO; m++) {
    	for (n = 0; n < N_ROW; n++) {
    	  prip_conc[m][n].new += prip_conc0[m][n].new;
    	}
  	  }

    for (m = 0; m < TYPY_PO; m++) {
      free((void *) prip_conc0[m]);
    }

    free((void **) prip_conc0);

 for (l = 0; l < N_ROW_PIC; l++) {                                                     /* preskalovani na potencial */
    celkovy_naboj = 0;
    i = l % Nx_UZLY;
    j = l / Nx_UZLY;
    real_cell_volume = Cell_volume;
    if((i==0)||(i==(Nx_UZLY-1))){
        if((j==0)||(j==Ny_UZLY_PIC-1)) {real_cell_volume = 0.25*Cell_volume;}
        else {real_cell_volume = 0.5*Cell_volume;}
    }
    if((j==0)||(j==(Ny_UZLY_PIC-1))){
        if((i==0)||(i==Nx_UZLY-1)) {real_cell_volume = 0.25*Cell_volume;}
        else {real_cell_volume = 0.5*Cell_volume;}
    }
    for (k = 0; k < TYPY_PO; k++) {
        celkovy_naboj += (naboj[k] * prip_conc[k][l].new);
        prip_conc[k][l].new /= real_cell_volume;
        prip_conc_PO[k][l].new   /= real_cell_volume;
     //prip_conc_porucha[k][l].new=prip_conc[k][l].new-conc_PO[k][l].new;
   }
   if((i==0)||(i==(Nx_UZLY-1))||(j==0)||(j==(Ny_UZLY_PIC-1))) celkovy_naboj=0;
   prip_b[l] = (-1/(eps0*Z_PO))*celkovy_naboj + prip_constpot[l];
 }

}


void VYPOCET_E(double *vyp_U, INTENZITA *vyp_E)          /* vypocet pole elektricke intenzity */
{ int i, j, k;
  int DOWN, UP;
																				/* vypocet Ex */
  UP = Nx_UZLY-1;

  for (j = 0; j < Ny_UZLY; j++) {
    k = j * Nx_UZLY;															/* leva hrana PO */
    vyp_E[k].x = - (vyp_U[k+1] - vyp_U[k]) / h_PO;
    for (i = 1; i < UP; i++) {													/* stred PO */
      k = j * Nx_UZLY + i;
      vyp_E[k].x = - (vyp_U[k+1] - vyp_U[k-1]) / (2 * h_PO);
    }

    k = (j + 1) * Nx_UZLY - 1;													/* prava hrana PO */
    vyp_E[k].x = - (vyp_U[k] - vyp_U[k-1]) / h_PO;
  }
																				/* vypocet Ey */
  for (k = 0; k < Nx_UZLY; k++) {
	vyp_E[k].y = - (vyp_U[k+Nx_UZLY] - vyp_U[k]) / h_PO;
  }

  UP = ((Ny_UZLY-1)*Nx_UZLY);

  for (k = Nx_UZLY; k < UP; k++) {												/* stred PO */
	vyp_E[k].y = - (vyp_U[k+Nx_UZLY] - vyp_U[k-Nx_UZLY]) / (2 * h_PO);
  }

  DOWN = ((Ny_UZLY-1)*Nx_UZLY);
  UP = (Nx_UZLY*Ny_UZLY);

  for (k = DOWN; k < UP; k++) {													/* horni hrana PO */
  	vyp_E[k].y = - (vyp_U[k] - vyp_U[k-Nx_UZLY]) / h_PO;
  }

}


void PRIPRAV_A(double *prip_Ax, int *prip_Ai, int *prip_Ap, int *prip_flag, double *prip_constpot)  /* pripravi pole Ax, Ai, Ap pro UMFPACK */
{                                                                               /* pole Ax - delka: pocet nenulovych prvku matice max = N_NOZERO */
 int i,j,k,l,m;                                                                   /* pole Ai - delka: pocet nenulovych prvku matice max = N_NOZERO */
 int p;                                /* aktualni pozice v poli Ax,Ai */       /* pole Ap - delka: pocet radku matice + 1,tj.N_ROW+1 */
 int r;                                /* aktualni pozice v poli Ap, tj.aktualni radek na kterem pocitam pocet nenulovych prvku */
 double x_uzel, y_uzel,x_r_help_pos,y_r_help_pos,x_r_help_neg,y_r_help_neg;
 double c,x_c,delta_x,y_c,delta_y;
 double pom_Ai[5], pom_Ax[5];													/* 0...spodni soused, 1...levy soused, 2...aktualni uzel, 3...pravy soused, 4...horni soused */
 double levy_soused, pravy_soused, horni_soused, dolni_soused;
 double radius = X_PO/2;

 prip_Ap[0] = 0;
 p = 0;

 for (k = 0; k < N_ROW_PIC; k++) {												/* cyklus pres vsechny body site tj. pres radky matice A */
   if (prip_flag[k] == 0) { prip_Ax[p] = 1.0 ;	prip_Ai[p] = k ; 	      p++;  /* aktualni bod site ma konstantni potencial */
							prip_Ap[k+1] = prip_Ap[k] + 1;
   }
   if (prip_flag[k] == 1) {
        i = k % Nx_UZLY;
        j = k / Nx_UZLY;
        /* uzel ma 4 sousedy, je uvnitr oblasti */
        if((i!=0)&&(i!=(Nx_UZLY-1))&&(j!=0)&&(j!=(Ny_UZLY-1))){
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - Nx_UZLY; p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - 1;       p++;
            prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + 1;       p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + Nx_UZLY; p++;
            prip_Ap[k+1] = prip_Ap[k] + 5;
        }
        else if (i==0){
                printf("nesmi nastat 1\n");
            if(j==0){ /* levy dolni roh */
                prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + 1;       p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + Nx_UZLY; p++;
                prip_Ap[k+1] = prip_Ap[k] + 3;
            }
            else if (j==(Ny_UZLY-1)){ /* levy horni roh*/
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - Nx_UZLY; p++;
                prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + 1;       p++;
                prip_Ap[k+1] = prip_Ap[k] + 3;
            }
            else { /* leva hrana */
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - Nx_UZLY; p++;
                prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + 1;       p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + Nx_UZLY; p++;
                prip_Ap[k+1] = prip_Ap[k] + 4;
            }

        }
        else if (i==(Nx_UZLY-1)){
                 printf("nesmi nastat 2\n");
            if(j==0){ /* pravy dolni roh */
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - 1;       p++;
                prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + Nx_UZLY; p++;
                prip_Ap[k+1] = prip_Ap[k] + 3;
            }
            else if (j==(Ny_UZLY-1)){ /* pravy horni roh */
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - Nx_UZLY; p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - 1;       p++;
                prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
                prip_Ap[k+1] = prip_Ap[k] + 3;
            }
            else { /* prava hrana */
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - Nx_UZLY; p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - 1;       p++;
                prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
                prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + Nx_UZLY; p++;
                prip_Ap[k+1] = prip_Ap[k] + 4;
            }
        }
        else if (j==0){ /* spodni hrana */
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - 1;       p++;
            prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + 1;       p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + Nx_UZLY; p++;
            prip_Ap[k+1] = prip_Ap[k] + 4;
        }
        else if (j==(Ny_UZLY-1)){ /* horni hrana */
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - Nx_UZLY; p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k - 1;       p++;
            prip_Ax[p] = -4.0;  prip_Ai[p] = k;           p++;
            prip_Ax[p] = 1.0 ;  prip_Ai[p] = k + 1;       p++;
            prip_Ap[k+1] = prip_Ap[k] + 4;
        }
   }

   //next for cycle will be commented
//   for(l=0;l<POCET_SOND;l++){
//    if (prip_flag[k] == l+2) { /* aktualni bod je u hranice l-te sondy, je vne */
//        /* urceni x,y souradnic uzlu */
//        r = 0;
//        i = k % Nx_UZLY;
//        j = k / Nx_UZLY;
//        x_uzel = X0_PO + i * h_PO;
//        y_uzel = Y0_PO + j * h_PO;
//        pom_Ax[2] = -4.0;  pom_Ai[2] = k;
//
//        /* ktery sousedni uzel lezi "v sonde"? */
//        levy_soused = sqrt(pow((x_uzel-h_PO)-X_S[l],2)+pow(y_uzel-Y_S[l],2));	/* vzdalenost leveho souseda od stredu sondy */
//        pravy_soused = sqrt(pow((x_uzel+h_PO)-X_S[l],2)+pow(y_uzel-Y_S[l],2));
//        horni_soused = sqrt(pow(x_uzel-X_S[l],2)+pow((y_uzel+h_PO)-Y_S[l],2));
//        dolni_soused = sqrt(pow(x_uzel-X_S[l],2)+pow((y_uzel-h_PO)-Y_S[l],2));
//
//        if (dolni_soused < R_S[l]) {
//            delta_y = sqrt(pow(R_S[l],2)-pow(X_S[l]-x_uzel,2));
//            y_c = Y_S[l] + delta_y;           	c = fabs(y_uzel-h_PO-y_c)/h_PO;
//            pom_Ax[4] = (1.0 - c / (2-c));  	pom_Ai[4] = k+Nx_UZLY;
//            pom_Ax[0] = 0;                  	pom_Ai[0] = 0;
//            prip_constpot[k] -= (2 * U_S[l]) / (2 - c);
//        }
//        else {
//            if (horni_soused < R_S[l]) {
//                delta_y = sqrt(pow(R_S[l],2)-pow(X_S[l]-x_uzel,2));
//                y_c = Y_S[l] - delta_y; 		c = fabs(y_uzel+h_PO-y_c)/h_PO;
//                pom_Ax[0] = (1.0 - c / (2-c));  pom_Ai[0] = k-Nx_UZLY;
//                pom_Ax[4] = 0; 					pom_Ai[4] = 0;
//                prip_constpot[k] -= (2 * U_S[l]) / (2 - c);
//            }
//            else {
//                pom_Ax[4] = 1.0;  				 	pom_Ai[4] = k + Nx_UZLY;
//                pom_Ax[0] = 1.0;                	pom_Ai[0] = k - Nx_UZLY;
//            }
//        }
//        if (levy_soused < R_S[l]) {
//            delta_x = sqrt(pow(R_S[l],2)-pow(Y_S[l]-y_uzel,2));
//            x_c = X_S[l] + delta_x; 		c = fabs(x_uzel-h_PO-x_c)/h_PO;
//            pom_Ax[3] = (1.0 - c / (2-c));  pom_Ai[3] = k+1;
//            pom_Ax[1] = 0; 					pom_Ai[1] = 0;
//            prip_constpot[k] -= (2 * U_S[l]) / (2 - c);
//        }
//        else {
//            if (pravy_soused < R_S[l]) {
//                delta_x = sqrt(pow(R_S[l],2)-pow(Y_S[l]-y_uzel,2));
//                x_c = X_S[l] - delta_x; 		c = fabs(x_uzel+h_PO-x_c)/h_PO;
//                pom_Ax[1] = (1.0 - c / (2-c)); 	pom_Ai[1] = k-1;
//                pom_Ax[3] = 0; 					pom_Ai[3] = 0;
//                prip_constpot[k] -= (2 * U_S[l]) / (2 - c);
//            }
//            else {
//                pom_Ax[1] = 1.0;  					pom_Ai[1] = k - 1;
//                pom_Ax[3] = 1.0;  					pom_Ai[3] = k + 1;
//            }
//        }
//        for (m = 0; m < 5; m++) {
//            if (pom_Ax[m] != 0) {
//                prip_Ax[p] = pom_Ax[m];
//                prip_Ai[p] = pom_Ai[m]; p++; r++;
//            }
//        }
//        prip_Ap[k+1] = prip_Ap[k] + r;
//    }
//   }
    // this part should be correct solution of the problem
    for(l=0;l<POCET_SOND;l++){
        if (prip_flag[k] == l+2) { /* aktualni bod je u hranice l-te sondy, je vne */
            /* urceni x,y souradnic uzlu */
            r = 0;
            i = k % Nx_UZLY;
            j = k / Nx_UZLY;
            // x and y coord of the k-th point
            x_uzel = X0_PO + i * h_PO;
            y_uzel = Y0_PO + j * h_PO;
            // y and x help value for calculation of the intersection with probe and since its sqrt we need both solutions
            x_r_help_pos = sqrt(pow(R_S[l],2)-pow((y_uzel-Y_S[l]),2)) + X_S[l];
            y_r_help_pos = sqrt(pow(R_S[l],2)-pow((x_uzel-X_S[l]),2)) + Y_S[l];
            x_r_help_neg = (-1)*sqrt(pow(R_S[l],2)-pow((y_uzel-Y_S[l]),2)) + X_S[l];
            y_r_help_neg = (-1)*sqrt(pow(R_S[l],2)-pow((x_uzel-X_S[l]),2)) + Y_S[l];
            // initialaze the lenghts to neighbours
            levy_soused  = 1.0;
            pravy_soused = 1.0;
            horni_soused = 1.0;
            dolni_soused = 1.0;
            //now the real values for lenghts !!but!! the change is made only if the nearest cell is inside the probe, first is left then right, down and up
            if(prip_flag[k-1] == 0)         levy_soused     = (x_uzel - x_r_help_pos) / h_PO;
            if(prip_flag[k+1] == 0)         pravy_soused    = (x_r_help_neg - x_uzel) / h_PO;
            if(prip_flag[k-Nx_UZLY] == 0)   dolni_soused    = (y_uzel - y_r_help_pos) / h_PO;
            if(prip_flag[k+Nx_UZLY] == 0)   horni_soused    = (y_r_help_neg - y_uzel) / h_PO;
            //and the final values
            prip_Ax[p] = 2.0 / (dolni_soused * (horni_soused+dolni_soused));    prip_Ai[p] = k - Nx_UZLY; p++;
            prip_Ax[p] = 2.0 / (levy_soused * (pravy_soused+levy_soused));      prip_Ai[p] = k - 1;       p++;
            prip_Ax[p] = (-2.0 * (pravy_soused*levy_soused+horni_soused*dolni_soused)) / (pravy_soused*levy_soused*horni_soused*dolni_soused);      prip_Ai[p] = k; p++;
            prip_Ax[p] = 2.0 / (pravy_soused * (pravy_soused+levy_soused));     prip_Ai[p] = k + 1;       p++;
            prip_Ax[p] = 2.0 / (horni_soused * (horni_soused+dolni_soused));    prip_Ai[p] = k + Nx_UZLY; p++;
            prip_Ap[k+1] = prip_Ap[k] + 5;
        }
    }
 }
}

 void BOUNDARY_GRADIENT(double *potential, double *pot_correction, double *edge_potential, double *charge_in_work_space)
 {
    // calculating Gauss law at the boundaries from gradient of potential + assumption that the front and back part will give 0 when summed
    double integral_left, integral_right, integral_down, integral_up, total_integral; //integrals on all sides of workspace
    double gradient, x_cord, y_cord;
    int i,j,k,l;
    double x,y,r;

    integral_left   = 0.0;
    integral_right  = 0.0;
    integral_down   = 0.0;
    integral_up     = 0.0;

    for(k=0; k<N_ROW; k++){
        pot_correction[k]=0.0;
    }

    //k=i+j*Nx_uzly
    for(j=0;j<Ny_UZLY;j++){
        //left side i=0
        k=j*Nx_UZLY;
        gradient = (potential[k+1]-potential[k])/h_PO;
        integral_left += gradient;
        //and right side i=Nx_uzly-1
        k=Nx_UZLY-1+j*Nx_UZLY;
        gradient = (potential[k]-potential[k-1])/h_PO;
        integral_right += gradient;
    }
    for(i=0;i<Nx_UZLY;i++){
        //down side j=0
        k=i;
        gradient = (potential[k+Nx_UZLY]-potential[k])/h_PO;
        integral_down += gradient;
        //and up side j=Ny_uzly-1
        k=i+(Ny_UZLY-1)*Nx_UZLY;
        gradient = (potential[k]-potential[k-Nx_UZLY])/h_PO;
        integral_up += gradient;
    }
    integral_left   *= Z_PO*Y_PO;
    integral_right  *= Z_PO*Y_PO;
    integral_down   *= Z_PO*X_PO;
    integral_up     *= Z_PO*X_PO;
    // need to calculate total_integral = integral through closed surface of -grad of potential => changed sign of the integral parts
    *charge_in_work_space = (integral_right-integral_left+integral_up-integral_down) / Real_Qpoint_value;

    //printf("naboj %e\n", *charge_in_work_space);

    for(i=0;i<Nx_UZLY;i++){
        //down side j=0
        k = i;
        l = i;
        x = X0_PO + i*h_PO - X_S[0];
        y = Y0_PO - Y_S[0];
        r = sqrt(x*x+y*y);
        pot_correction[l] = *charge_in_work_space * (h_PO/(r+h_PO)) * edge_potential[k];
        //and up side j=Ny_uzly-1
        k=i+Nx_UZLY+Ny_UZLY-2;
        l=i+(Ny_UZLY-1)*Nx_UZLY;
        y = Y0_PO + (Ny_UZLY-1)*h_PO - Y_S[0];
        r = sqrt(x*x+y*y);
        pot_correction[l] = *charge_in_work_space * (h_PO/(r+h_PO)) * edge_potential[k];
    }
    for(j=1;j<Ny_UZLY-1;j++){
        //left side i=0
        k=j-1+Nx_UZLY;
        l=j*Nx_UZLY;
        x = X0_PO - X_S[0];
        y = Y0_PO + j*h_PO - Y_S[0];
        r = sqrt(x*x+y*y);
        pot_correction[l] = *charge_in_work_space * (h_PO/(r+h_PO)) * edge_potential[k];
        //and right side i=Nx_uzly-1
        k=j+2*Nx_UZLY+Ny_UZLY-3;
        l=Nx_UZLY-1+j*Nx_UZLY;
        x = X0_PO + (Nx_UZLY-1)*h_PO - X_S[0];
        r = sqrt(x*x+y*y);
        pot_correction[l] = *charge_in_work_space * (h_PO/(r+h_PO)) * edge_potential[k];
    }
 }

// function to crate field of potential at the boundary of workspace for point charge at the probe position without Q/eps0 in the eq
// number of nodes in edge potential is N_BC_NODES and it should start at left bottom corner of workspace and go counter-clockwise
// data stored for bottom line from left to right - Nx_Uzly times, for right side from down to up - Ny_UZLY-2 times
// upper side from left to right again - Nx_Uzly times, for fest side again from down to up - Ny_UZLY-2 times
 void POINT_CHARGE_EDGE_POTENTIAL(double *edge_potential){
    int k;
    double x_from_probe, y_from_probe;

    //bottom
    for(k=0; k<Nx_UZLY;k++){
        //i=k; j=0;
        x_from_probe = X0_PO + k * h_PO - X_S[0];
        y_from_probe = Y0_PO - Y_S[0];
        edge_potential[k] = POTENTIAL_OF_POINT_CHARGE(&x_from_probe,&y_from_probe);
    }
    //right
    for(k=Nx_UZLY; k<Nx_UZLY+Ny_UZLY-2;k++){
        //i=Nx_UZLY-1; j=k-Nx_UZLY+1;
        x_from_probe = X0_PO + (Nx_UZLY-1) * h_PO - X_S[0];
        y_from_probe = Y0_PO + (k-Nx_UZLY+1) * h_PO - Y_S[0];
        edge_potential[k] = POTENTIAL_OF_POINT_CHARGE(&x_from_probe,&y_from_probe);
    }
    //upper
    for(k=Nx_UZLY+Ny_UZLY-2; k<2*Nx_UZLY+Ny_UZLY-2;k++){
        //i=k-(Nx_UZLY+Ny_UZLY-2); j=Ny_UZLY-1;
        x_from_probe = X0_PO + (k-(Nx_UZLY+Ny_UZLY-2)) * h_PO - X_S[0];
        y_from_probe = Y0_PO + (Ny_UZLY-1) * h_PO - Y_S[0];
        edge_potential[k] = POTENTIAL_OF_POINT_CHARGE(&x_from_probe,&y_from_probe);
    }
    //left
    for(k=2*Nx_UZLY+Ny_UZLY-2; k<2*Nx_UZLY+2*Ny_UZLY-4; k++){
        //i=0; j=k-(2*Nx_UZLY+Ny_UZLY-2)+1;
        x_from_probe = X0_PO - X_S[0];
        y_from_probe = Y0_PO + (k-(2*Nx_UZLY+Ny_UZLY-2)+1) * h_PO - Y_S[0];
        edge_potential[k] = POTENTIAL_OF_POINT_CHARGE(&x_from_probe,&y_from_probe);
    }
 }
//point charge where Q=1
double POTENTIAL_OF_POINT_CHARGE(double *x, double *y){
    double r;
    r = sqrt(( *x - X_S[0] )*( *x - X_S[0] )+( *y - Y_S[0] )*( *y - Y_S[0] ));
    return 1/(4*M_PI*eps0*r);
}

double Q_PART_OF_POINTCHARGE(){
    int i,j;
    double x,y,grad_left,grad_right,grad_down,grad_up,integral,delta;

    grad_left = 0.0;
    grad_right = 0.0;
    grad_down = 0.0;
    grad_up = 0.0;
    integral = 0.0;

    for(j=0;j<Ny_UZLY;j++){
        y = Y0_PO + j*h_PO - Y_S[0] ;
    //left
        x = X0_PO - X_S[0];
        delta = x+h_PO;
        grad_left += (POTENTIAL_OF_POINT_CHARGE(&delta,&y)-POTENTIAL_OF_POINT_CHARGE(&x,&y))/h_PO;
    //right
        x = X0_PO + X_PO - X_S[0];
        delta = x-h_PO;
        grad_right += (POTENTIAL_OF_POINT_CHARGE(&x,&y)-POTENTIAL_OF_POINT_CHARGE(&delta,&y))/h_PO;
    }
    for(i=0;i<Nx_UZLY;i++){
        x = X0_PO +i*h_PO - X_S[0];
    //down
        y = Y0_PO - Y_S[0];
        delta = y+h_PO;
        grad_down += (POTENTIAL_OF_POINT_CHARGE(&x,&delta)-POTENTIAL_OF_POINT_CHARGE(&x,&y))/h_PO;
    //up
        y = Y0_PO + Y_PO - Y_S[0];
        delta = y-h_PO;
        grad_up += (POTENTIAL_OF_POINT_CHARGE(&x,&y)-POTENTIAL_OF_POINT_CHARGE(&x,&delta))/h_PO;
    }

    grad_left *= Z_PO*Y_PO;
    grad_right *= Z_PO*Y_PO;
    grad_down *= Z_PO*X_PO;
    grad_up *= Z_PO*X_PO;
    integral = (grad_right-grad_left+grad_up-grad_down); // integral is Q_part from point charge of Q=1 (here i dont need *eps0 because i dont use it upper neither)
    return integral;
}
