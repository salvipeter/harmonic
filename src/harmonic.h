#ifndef HARMONIC_H
#define HARMONIC_H

/* Based on: T. DeRose, M. Meyer: Harmonic Coordinates. Pixar Animation Studios */

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct HarmonicMap *harmonic_create(const double *min, const double *max, size_t levels);

void harmonic_add_point(struct HarmonicMap *map, const double *point);

void harmonic_add_line(struct HarmonicMap *map, const double *from, const double *to);

void harmonic_add_curve(struct HarmonicMap *map, const double *points, size_t n, size_t resolution);

void harmonic_solve(struct HarmonicMap *map, double epsilon, bool biharmonic);

bool harmonic_eval(const struct HarmonicMap *map, const double *point, double *value);

void harmonic_write_ppm(const struct HarmonicMap *map, const char *filename);

void harmonic_free(struct HarmonicMap *map);

#ifdef __cplusplus
}
#endif

#endif
