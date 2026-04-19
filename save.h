#ifndef SAVE_H_INCLUDED
#define SAVE_H_INCLUDED

extern void save_n_E_U(double *a_U_old, INTENZITA *a_E_old, KONCENTRACE **a_conc, KONCENTRACE **a_conc_PO, KONCENTRACE **a_conc_porucha, int a_iterace, int kdy);
extern void save_param(DDCASTICE *save_oblast, DDCASTICE *save_zdroj_1, DDCASTICE *save_zdroj_2, DDCASTICE *save_zdroj_3, DDCASTICE *save_zdroj_4, int N_castic_PO, int N_castic_zdroj, int aktualni_iterace);
extern void save_param_1(DDCASTICE *save_oblast, DDCASTICE *save_zdroj_1, DDCASTICE *save_zdroj_2, DDCASTICE *save_zdroj_3, DDCASTICE *save_zdroj_4, int N_castic_PO, int N_castic_zdroj, int aktualni_iterace);
extern void save_model_parameters(FILE *mp_file);

#endif
