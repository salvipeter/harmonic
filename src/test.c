#include "harmonic.h"

int main(int argc, char **argv) {
  double input[] = {0, 0, 0,
                    6, 0, 0,
                    6, 6, 0,
                    4, 6, 0,
                    4, 2, 0,
                    2, 2, 0,
                    2, 6, 1,
                    0, 6, 1};
  struct HarmonicMap *map = harmonic_init(8, input, 9, 1.0e-5);
  double point[] = { 0.5, 0.5 };
  harmonic_eval(map, point);
  harmonic_free(map);
  return 0;
}
