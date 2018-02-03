#ifndef HARMONIC_H
#define HARMONIC_H

/* Based on: T. DeRose, M. Meyer: Harmonic Coordinates. Pixar Animation Studios */

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* points contains 3*size values: x0,y0,v0, x1,y1,v1, etc.
   epsilon is the minimal average change to continue the solver iteration */
struct HarmonicMap *harmonic_init(unsigned int size, double *points, unsigned int levels,
                                  double epsilon);

bool harmonic_eval(struct HarmonicMap *map, double *point, double *value);

void harmonic_write_ppm(struct HarmonicMap *map, char *filename);

void harmonic_free(struct HarmonicMap *map);

/* Returns the number of vertices */
unsigned int harmonic_mesh_size(struct HarmonicMap *map, unsigned int downsampling);

/* vertices should be 2*mesh_size long, triangles 6*mesh_size
   the return value is the real number of triangles */
unsigned int harmonic_mesh(struct HarmonicMap *map, unsigned int downsampling,
                           double *vertices, unsigned int *triangles);

#ifdef __cplusplus
}
#endif

#endif
