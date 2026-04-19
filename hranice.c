#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "global.h"
#include "hranice.h"

int HRANICNI_PODMINKY(DDCASTICE *castice)
{
    int i;
     /* kontrola opusteni pracovni oblasti, pokud ne priradi se vypadla=0 */

    if ( (castice->x[0]>(-(X_PO/2.0))) && (castice->x[0]<(X_PO/2.0)) && (castice->x[1]>(-(Y_PO/2.0))) && (castice->x[1]<(Y_PO/2.0)) )
    {
        for(i=0; i<POCET_SOND; i++){
            if ((sqrt(pow(castice->x[0]-X_S[i],2)+pow(castice->x[1]-Y_S[i],2))<=R_S[i])&&(SONDA[i]==1))  return (i+2);
        }
        return 0;
    }
    else                                                                                                 return 1;
}
