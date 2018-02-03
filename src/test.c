#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

int main(int argc, char **argv) {
  double input[] = { 0, 0, 0,
                     6, 0, 0,
                     6, 6, 0,
                     3.5, 6, 0,
                     3, 3, 1,
                     2.5, 6, 0,
                     0, 6, 0 };
  bool success;
  double result, point[] = { 2, 2 };
  unsigned int levels = 9;
  struct HarmonicMap *map;

  map = harmonic_init(7, input, levels, 1.0e-5);

  /* Evaluation test */
  success = harmonic_eval(map, point, &result);
  if (success)
    printf("%lf\n", result);
  else
    printf("N/A\n");

  /* PPM output test */
  harmonic_write_ppm(map, "/tmp/test.ppm");

  /* Mesh generation test */

  unsigned int v_size = harmonic_mesh_size(map, 2);
  double *v = (double *)malloc(sizeof(double) * v_size * 2);
  unsigned int *t = (unsigned int *)malloc(sizeof(unsigned int) * v_size * 6);
  unsigned int t_size = harmonic_mesh(map, 2, v, t);

  FILE *f = fopen("/tmp/test.obj", "w");
  for (unsigned int i = 0; i < v_size; ++i)
    fprintf(f, "v %f %f 0.0\n", v[2*i], v[2*i+1]);
  for (unsigned int i = 0; i < t_size; ++i)
    fprintf(f, "f %d %d %d\n", t[3*i] + 1, t[3*i+1] + 1, t[3*i+2] + 1);
  fclose(f);

  free(t);
  free(v);

  harmonic_free(map);

  return 0;
}
