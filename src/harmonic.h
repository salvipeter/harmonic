#ifndef HARMONIC_H
#define HARMONIC_H

/* Based on: T. DeRose, M. Meyer: Harmonic Coordinates. Pixar Animation Studios */

struct GridValue {
  enum GridValueType { UNTYPED, EXTERIOR, BOUNDARY, INTERIOR } type;
  double value;
};

struct HarmonicMap {
  unsigned int size;
  struct GridValue *grid;
  double offset[2];
  double scaling;
};

/* points contains 3*size values: x0,y0,v0, x1,y1,v1, etc.
   epsilon is the minimal average change to continue the solver iteration */
struct HarmonicMap *harmonic_init(unsigned int size, double *points, unsigned int levels,
                                  double epsilon);

double harmonic_eval(struct HarmonicMap *map, double *point);

void harmonic_free(struct HarmonicMap *map);

#endif
