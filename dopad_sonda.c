#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

#include "global.h"
#include "dopad_sonda.h"

int SGN(double argument){
    int vystup;

    if      (argument>0)    vystup=1;
    else if (argument<0)    vystup=(-1);
    else                    vystup=0;

    return vystup;
}

int DOPAD_SONDA(DDCASTICE *castice, int cislo_sondy){

    double  a,b,c,diskriminant,ksi_help_1,ksi_help_2,ksi,x_help,y_help,x_kruh,y_kruh,fi,fi_deg;
    int     vysledek;

    x_help=castice->x[0]-X_S[cislo_sondy];
    y_help=castice->x[1]-Y_S[cislo_sondy];

    a=pow(castice->v[0],2)+pow(castice->v[1],2);
    b=(-2)*(x_help*castice->v[0]+y_help*castice->v[1]);
    c=x_help*x_help+y_help*y_help-R_S[cislo_sondy]*R_S[cislo_sondy];

    diskriminant=b*b-4*a*c;
    if(diskriminant<0)("printf error diskriminant zaporny \n");

    ksi_help_1=((-b)+sqrt(diskriminant))/(2*a);
    ksi_help_2=((-b)-sqrt(diskriminant))/(2*a);

    if(ksi_help_1<=DT[castice->typ])    ksi=ksi_help_1;
    else                                ksi=ksi_help_2;

    x_kruh=(castice->x[0]-castice->v[0]*ksi)-X_S[cislo_sondy];
    y_kruh=(castice->x[1]-castice->v[1]*ksi)-Y_S[cislo_sondy];

    fi = atan2(y_kruh,x_kruh);
    if (fi<0.0) fi+=2*M_PI;

    fi_deg = fi*180/M_PI;
    vysledek = (int) round(fi_deg);
    if(vysledek==360) vysledek=0;

    return vysledek;
}
