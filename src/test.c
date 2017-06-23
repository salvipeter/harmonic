#include "harmonic.h"

int main(int argc, char **argv) {
  double input[] = { 2, 0, 1,
                     4, 1, 1,
                     6, 0, 0,
                     7, 2, 0,
                     5, 3, 0,
                     6, 4, 0,
                     7, 7, 0,
                     5, 5, 0,
                     4, 7, 0,
                     2, 5, 0,
                     0, 7, 0,
                     1, 4, 0,
                     3, 3, 0,
                     1, 2, 0};
  struct HarmonicMap *map = harmonic_init(14, input, 9, 1.0e-5);
  double point[] = { 0.5, 0.5 };
  harmonic_eval(map, point);
  harmonic_write_ppm(map, "/tmp/csillag.ppm");
  harmonic_free(map);
  return 0;
}
