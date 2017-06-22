#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN_LEVEL 3

void flood_fill(struct HarmonicMap *map, unsigned int n, unsigned int x, unsigned int y) {
  if (map->grid[y*n+x].type != UNTYPED)
    return;
  map->grid[y*n+x].type = EXTERIOR;
  if (x > 0)
    flood_fill(map, n, x - 1, y);
  if (x < n - 1)
    flood_fill(map, n, x + 1, y);
  if (y > 0)
    flood_fill(map, n, x, y - 1);
  if (y < n - 1)
    flood_fill(map, n, x, y + 1);
}

void solve(struct GridValue *grid, unsigned int level, double epsilon) {
  unsigned int n = pow(2, level);
  if (level > MIN_LEVEL) {
    /* Generate a coarser grid and solve that first to get good starting values */
    unsigned int level1 = level - 1, n1 = pow(2, level1);
    struct GridValue *grid1 = (struct GridValue *)malloc(n1 * n1 * sizeof(struct GridValue));
    for (unsigned int i = 0; i < n1; ++i)
      for (unsigned int j = 0; j < n1; ++j) {
        if (grid[2*j*n+2*i].type == BOUNDARY ||
            (i < n1 && grid[2*j*n+2*i+1].type == BOUNDARY) ||
            (j < n1 && grid[(2*j+1)*n+2*i].type == BOUNDARY) ||
            (i < n1 && j < n1 && grid[(2*j+1)*n+2*i+1].type == BOUNDARY)) {
          grid1[j*n1+i].type = BOUNDARY;
          grid1[j*n1+i].value = 0.0;
          int count = 0;
          if (grid[2*j*n+2*i].type == BOUNDARY) {
            ++count;
            grid1[j*n1+i].value += grid[2*j*n+2*i].value;
          }
          if (i < n1 && grid[2*j*n+2*i+1].type == BOUNDARY) {
            ++count;
            grid1[j*n1+i].value += grid[2*j*n+2*i+1].value;
          }
          if (j < n1 && grid[(2*j+1)*n+2*i].type == BOUNDARY) {
            ++count;
            grid1[j*n1+i].value += grid[(2*j+1)*n+2*i].value;
          }
          if (i < n1 && j < n1 && grid[(2*j+1)*n+2*i+1].type == BOUNDARY) {
            ++count;
            grid1[j*n1+i].value += grid[(2*j+1)*n+2*i+1].value;
          }
          if (count > 0)
            grid1[j*n1+i].value /= (double)count;
        } else if (grid[2*j*n+2*i].type == EXTERIOR)
          grid1[j*n1+i].type = EXTERIOR;
        else {
          grid1[j*n1+i].type = INTERIOR;
          grid1[j*n1+i].value = 0.0;
        }
      }
    solve(grid1, level1, epsilon);
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j)
        if (grid[j*n+i].type == INTERIOR)
          grid[j*n+i].value = grid1[(j/2)*n1+i/2].value;
    free(grid1);
  }

  /* Solve by iteration */
  double change;
  unsigned int count, iteration = 0;
  struct GridValue *tmp = grid;
  for (unsigned int i = 0; i < n * n; ++i)
    tmp[i].type = grid[i].type;
  do {
    ++iteration;
    change = 0.0;
    count = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j)
        if (grid[j*n+i].type == INTERIOR) {
          int neighbors = 0;
          tmp[j*n+i].value = 0.0;
          if (j > 0 && grid[(j-1)*n+i].type != EXTERIOR) {
            ++neighbors;
            tmp[j*n+i].value += grid[(j-1)*n+i].value;
          }
          if (i > 0 && grid[j*n+i-1].type != EXTERIOR) {
            ++neighbors;
            tmp[j*n+i].value += grid[j*n+i-1].value;
          }
          if (j < n && grid[(j+1)*n+i].type != EXTERIOR) {
            ++neighbors;
            tmp[j*n+i].value += grid[(j+1)*n+i].value;
          }
          if (i < n && grid[j*n+i+1].type != EXTERIOR) {
            ++neighbors;
            tmp[j*n+i].value += grid[j*n+i+1].value;
          }
          assert(neighbors > 0);
          tmp[j*n+i].value /= (double)neighbors;
          ++count;
          change += fabs(tmp[j*n+i].value - grid[j*n+i].value);
        }
    change /= (double)count;
  } while (change > epsilon);

  /* Write a PPM file for testing */
  char str[255];
  sprintf(str, "/tmp/proba-%d.ppm", level);
  FILE *f = fopen(str, "w");
  fprintf(f, "P3\n%d %d\n255\n", n, n);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j)
      if (grid[j*n+i].type == EXTERIOR)
        fprintf(f, "255 0 0 ");
      else
        fprintf(f, "0 0 %d ", (int)(grid[j*n+i].value * 255.0));
    fprintf(f, "\n");
  }
  fclose(f);
}

struct HarmonicMap *harmonic_init(unsigned int size, double *points, unsigned int levels,
                                  double epsilon) {
  /* Create the grid */
  unsigned int n = pow(2, levels);
  struct HarmonicMap *map = (struct HarmonicMap *)malloc(sizeof(struct HarmonicMap));
  map->size = n;
  map->grid = (struct GridValue *)malloc(n * n * sizeof(struct GridValue));

  /* Bounding box computation */
  double min[2], max[2];
  min[0] = max[0] = points[0]; min[1] = max[1] = points[1];
  for (unsigned int i = 1; i < size; ++i) {
    if (points[i*3] < min[0])
      min[0] = points[i*3];
    else if (points[i*3] > max[0])
      max[0] = points[i*3];
    if (points[i*3+1] < min[1])
      min[1] = points[i*3+1];
    else if (points[i*3+1] > max[1])
      max[1] = points[i*3+1];
  }
  double length = MAX(max[0] - min[0], max[1] - min[1]);

  /* Add a margin of 2.5% on all sides */
  map->offset[0] = min[0] - length * 0.025;
  map->offset[1] = min[1] - length * 0.025;
  map->scaling = length * 1.05 / (double)n;

  /* Fill boundary cells */
  for (unsigned int i = 0; i < n * n; ++i)
    map->grid[i].type = UNTYPED;
  int x0 = (points[size*3-3] - map->offset[0]) / map->scaling;
  int y0 = (points[size*3-2] - map->offset[1]) / map->scaling;
  double v0 = points[size*3-1];
  for (unsigned int i = 0; i < size; ++i) {
    /* line drawing from Rosetta Code [Bresenham's algorithm] */
    int x1 = (points[i*3] - map->offset[0]) / map->scaling;
    int y1 = (points[i*3+1] - map->offset[1]) / map->scaling;
    double v1 = points[i*3+2];
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1; 
    int err = (dx > dy ? dx : -dy) / 2, e2;
    for(;;) {
      double ratio;             /* linear interpolation along the sides */
      if (err > 0)
        ratio = (double)abs(x1 - x0) / (double)dx;
      else
        ratio = (double)abs(y1 - y0) / (double)dy;
      map->grid[y0*n+x0].type = BOUNDARY;
      map->grid[y0*n+x0].value = v0 * ratio + v1 * (1.0 - ratio);
      if (x0 == x1 && y0 == y1) break;
      e2 = err;
      if (e2 > -dx) { err -= dy; x0 += sx; }
      if (e2 <  dy) { err += dx; y0 += sy; }
    }
    v0 = v1;
  }

  /* Fill exterior and interior */
  flood_fill(map, n, 0, 0);
  for (unsigned int i = 0; i < n * n; ++i)
    if (map->grid[i].type == UNTYPED) {
      map->grid[i].type = INTERIOR;
      map->grid[i].value = 0.0;
    }

  /* Find the solution for the discrete problem */
  solve(map->grid, levels, epsilon);

  return map;
}

double harmonic_eval(struct HarmonicMap *map, double *point) {
  double x = (point[0] - map->offset[0]) / map->scaling;
  double y = (point[1] - map->offset[1]) / map->scaling;
  int i = (int)x, j = (int)y;
  if (i == map->size - 1)
    i = map->size - 2;
  if (j == map->size - 1)
    j = map->size - 2;
  return map->grid[j*n+i] * (1.0 - y + j) * (1.0 - x + i) +
    map->grid[(j+1)*n+i] * (y - j) * (1.0 - x + i) +
    map->grid[j*n+i+1] * (1.0 - y + j) * (x - i) +
    map->grid[(j+1)*n+i+1] * (y - j) * (x - i);
}

void harmonic_free(struct HarmonicMap *map) {
  free(map->grid);
  free(map);
}
