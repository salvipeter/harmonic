#ifndef HARMONIC_H
#define HARMONIC_H

/* Based on: T. DeRose, M. Meyer: Harmonic Coordinates. Pixar Animation Studios */

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* points contains 3*size values: x0,y0,v0, x1,y1,v1, etc.
   epsilon is the minimal average change to continue the solver iteration */
struct HarmonicMap *harmonic_init(size_t size, double *points, size_t levels, double epsilon,
                                  bool biharmonic);

bool harmonic_eval(struct HarmonicMap *map, double *point, double *value);

void harmonic_write_ppm(struct HarmonicMap *map, char *filename);

void harmonic_free(struct HarmonicMap *map);

#ifdef __cplusplus
}
#endif

#endif
