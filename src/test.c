#include <math.h>
#include <stdio.h>

#include "harmonic.h"

int main(int argc, char **argv) {
  double input[] = { 0, 0, 0,
                     6, 0, 0,
                     6, 6, 0,
                     4, 6, 0,
                     4, 2, 0,
                     2, 2, 0,
                     2, 6, 1,
                     0, 6, 0 };
  unsigned int levels = 9, n = pow(2, levels);
  struct HarmonicMap *map, *map2;
  FILE *f;

  map = harmonic_init(8, input, levels, 1.0e-5);
  input[6*3+2] = 0; input[7*3+2] = 1;
  map2 = harmonic_init(8, input, levels, 1.0e-5);

  f = fopen("/tmp/test.ppm", "w");
  fprintf(f, "P3\n%d %d\n255\n", n, n);
  double x, y;
  for (unsigned int i = 0; i < n; ++i) {
    x = 6.0 * (double)i / n;
    for (unsigned int j = 0; j < n; ++j) {
      y = 6.0 * (double)j / n;
      double p[] = {x, y}, a, b;
      if (harmonic_eval(map, p, &a) && harmonic_eval(map2, p, &b)) {
        if (a + b > 1.0e-5) {
          double alpha = a / (a + b);
          fprintf(f, "0 0 %d ", (int)(alpha * 255.0));
        } else {
          fprintf(f, "0 255 0 ");
        }
      } else
        fprintf(f, "255 0 0 ");
    }
    fprintf(f, "\n");
  }
  fclose(f);

  harmonic_free(map2);
  harmonic_free(map);
  return 0;
}
