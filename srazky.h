#ifndef SRAZKY_H_INCLUDED
#define SRAZKY_H_INCLUDED

extern void srazkove_procesy(DDCASTICE *castice, int thread_n);
extern double scatter_electron_majorant();
extern void open_cross_sections_file(int *counter_el, int *counter_ex, int *counter_ion);
extern void read_el_cross_sections(double **sigma_el, double **sigma_ex, double **sigma_ion);
extern double cross_section_interpolate(double **sigma_array, int array_size, double particle_energy);

#endif // SRAZKY_H_INCLUDED
