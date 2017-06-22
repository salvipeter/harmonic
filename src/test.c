#include "harmonic.h"

int main(int argc, char **argv) {
  double input[] = {0, 0, 0,
                    5, 0, 0,
                    3, 3, 1,
                    5, 5, 1,
                    0, 5, 0,
                    1, 2, 0};
  struct HarmonicMap *map = harmonic_init(8, input, 9, 1.0e-5);
  double point[] = { 0.5, 0.5 };
  harmonic_eval(map, point);
  harmonic_free(map);
  return 0;
}
