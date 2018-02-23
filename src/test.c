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

  map = harmonic_init(7, input, levels, 1.0e-5, true);

  /* Evaluation test */
  success = harmonic_eval(map, point, &result);
  if (success)
    printf("%lf\n", result);
  else
    printf("N/A\n");

  /* PPM output test */
  harmonic_write_ppm(map, "/tmp/test.ppm");

  harmonic_free(map);

  return 0;
}
