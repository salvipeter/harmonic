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
  size_t levels = 9;
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

  size_t v_size = harmonic_mesh_size(map, 2);
  double *v = (double *)malloc(sizeof(double) * v_size * 2);
  size_t *t = (size_t *)malloc(sizeof(size_t) * v_size * 6);
  size_t t_size = harmonic_mesh(map, 2, v, t);

  FILE *f = fopen("/tmp/test.obj", "w");
  for (size_t i = 0; i < v_size; ++i)
    fprintf(f, "v %f %f 0.0\n", v[2*i], v[2*i+1]);
  for (size_t i = 0; i < t_size; ++i)
    fprintf(f, "f %zu %zu %zu\n", t[3*i] + 1, t[3*i+1] + 1, t[3*i+2] + 1);
  fclose(f);

  free(t);
  free(v);

  harmonic_free(map);

  return 0;
}
