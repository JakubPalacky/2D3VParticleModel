#ifndef PIC_H
#define PIC_H

extern void FLAG_CONSTPOT(int *prip_flag, double *prip_constpot);
extern void CONSTPOT_BOUNDARY_CHANGE(double *prip_constpot_boundary, double *pole_constpot, double *pole_U_BC);
extern void PRIPRAV_B(DDCASTICE *prip_PO, DDCASTICE *prip_porucha, double *prip_b, double *prip_constpot, KONCENTRACE **prip_conc,KONCENTRACE **prip_conc_PO, KONCENTRACE **prip_conc_porucha);
extern void VYPOCET_E(double *vyp_U, INTENZITA *vyp_E);
extern void PRIPRAV_A(double *prip_Ax, int *prip_Ai, int *prip_Ap, int *prip_flag, double *prip_constpot);
extern void BOUNDARY_GRADIENT(double *potential, double *pot_correction, double *edge_potential, double *charge_in_work_space);
extern void POINT_CHARGE_EDGE_POTENTIAL(double *edge_potential);
extern double POTENTIAL_OF_POINT_CHARGE(double *x, double *y);
extern double Q_PART_OF_POINTCHARGE();

#endif
