#include <math.h>
#include <stdio.h>

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
  success = harmonic_eval(map, point, &result);
  if (success)
    printf("%lf\n", result);
  else
    printf("N/A\n");
  harmonic_write_ppm(map, "/tmp/test.ppm");
  harmonic_free(map);
  return 0;
}
