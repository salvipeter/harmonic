#ifndef HARMONIC_H
#define HARMONIC_H

/* Based on: T. DeRose, M. Meyer: Harmonic Coordinates. Pixar Animation Studios */

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct HarmonicMap *harmonic_create(double *min, double *max, size_t levels);

void harmonic_add_point(struct HarmonicMap *map, double *point);

void harmonic_add_line(struct HarmonicMap *map, double *from, double *to);

void harmonic_solve(struct HarmonicMap *map, double epsilon, bool biharmonic);

bool harmonic_eval(struct HarmonicMap *map, double *point, double *value);

void harmonic_write_ppm(struct HarmonicMap *map, char *filename);

void harmonic_free(struct HarmonicMap *map);

#ifdef __cplusplus
}
#endif

#endif
