#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

struct GridValue {
  bool boundary;
  double value;
};

struct HarmonicMap {
  size_t levels, size;
  struct GridValue *grid;
  double offset[2];
  double scaling;
};

#define MIN_LEVEL 3

#define MAX(a,b) (((a)>(b))?(a):(b))

void solveHarmonic(struct GridValue *grid, size_t n, double epsilon) {
  double change;
  do {
    change = 0.0;
    size_t count = 0, index = n + 1;
    for (size_t j = 1, n_1 = n - 1; j < n_1; ++j) {
      for (size_t i = 1; i < n_1; ++i, ++index)
        if (!grid[index].boundary) {
          double value = 0.0;
          value += grid[index-n].value;
          value += grid[index-1].value;
          value += grid[index+n].value;
          value += grid[index+1].value;
          value /= 4.0;
          change += fabs(grid[index].value - value);
          grid[index].value = value;
          ++count;
        }
      index += 2;
    }
    /* Boundary cases: not handled [the result is the same] */
    change /= (double)count;
  } while (change > epsilon);
}

void solveBiharmonic(struct GridValue *grid, size_t n, double epsilon) {
  double change;
  do {
    change = 0.0;
    size_t count = 0, n_2 = n - 2, n2 = n * 2, index = n2 + 2;
    for (size_t j = 2; j < n_2; ++j) {
      for (size_t i = 2; i < n_2; ++i, ++index)
        if (!grid[index].boundary) {
          double value = 0.0;
          value += grid[index-n].value;
          value += grid[index-1].value;
          value += grid[index+n].value;
          value += grid[index+1].value;
          value *= 4.0;
          value -= grid[index-n-1].value;
          value -= grid[index-n+1].value;
          value -= grid[index+n-1].value;
          value -= grid[index+n+1].value;
          value *= 2.0;
          value -= grid[index-n2].value;
          value -= grid[index-2].value;
          value -= grid[index+n2].value;
          value -= grid[index+2].value;
          value /= 20.0;
          change += fabs(grid[index].value - value);
          grid[index].value = value;
          ++count;
        }
      index += 4;
    }
    /* Boundary cases */
    for (size_t j = 0; j < n; ++j)
      for (size_t i = 0; i < n; ++i) {
        if (j >= 2 && j < n_2 && i >= 2 && i < n_2)
          continue;
        size_t index = j * n + i;
        if (grid[index].boundary)
          continue;
        double weight = 0.0;
        double value = 0.0;
        if (i > 0) {
          size_t k = index - 1;
          weight += 2.0;
          value += grid[k].value * 2.0;
          int neighbors = 4;
          if (i == 1)
            --neighbors;
          if (j == 0)
            --neighbors;
          if (j == n - 1)
            --neighbors;
          double w = 1.0 / (double)neighbors;
          if (i > 1) {
            weight -= w;
            value -= grid[k-1].value * w;
          }
          if (j > 0) {
            weight -= w;
            value -= grid[k-n].value * w;
          }
          if (j < n - 1) {
            weight -= w;
            value -= grid[k+n].value * w;
          }
        }
        if (i < n - 1) {
          size_t k = index + 1;
          weight += 2.0;
          value += grid[k].value * 2.0;
          int neighbors = 4;
          if (i == n - 2)
            --neighbors;
          if (j == 0)
            --neighbors;
          if (j == n - 1)
            --neighbors;
          double w = 1.0 / (double)neighbors;
          if (i < n - 2) {
            weight -= w;
            value -= grid[k+1].value * w;
          }
          if (j > 0) {
            weight -= w;
            value -= grid[k-n].value * w;
          }
          if (j < n - 1) {
            weight -= w;
            value -= grid[k+n].value * w;
          }
        }
        if (j > 0) {
          size_t k = index - n;
          weight += 2.0;
          value += grid[k].value * 2.0;
          int neighbors = 4;
          if (i == 0)
            --neighbors;
          if (i == n - 1)
            --neighbors;
          if (j == 1)
            --neighbors;
          double w = 1.0 / (double)neighbors;
          if (i > 0) {
            weight -= w;
            value -= grid[k-1].value * w;
          }
          if (i < n - 1) {
            weight -= w;
            value -= grid[k+1].value * w;
          }
          if (j > 1) {
            weight -= w;
            value -= grid[k-n].value * w;
          }
        }
        if (j < n - 1) {
          size_t k = index + n;
          weight += 2.0;
          value += grid[k].value * 2.0;
          int neighbors = 4;
          if (i == 0)
            --neighbors;
          if (i == n - 1)
            --neighbors;
          if (j == n - 2)
            --neighbors;
          double w = 1.0 / (double)neighbors;
          if (i > 0) {
            weight -= w;
            value -= grid[k-1].value * w;
          }
          if (i < n - 1) {
            weight -= w;
            value -= grid[k+1].value * w;
          }
          if (j < n - 2) {
            weight -= w;
            value -= grid[k+n].value * w;
          }
        }
        value /= weight;
        change += fabs(grid[index].value - value);
        grid[index].value = value;
        ++count;
      }
    change /= (double)count;
  } while (change > epsilon);
}

void solve(struct GridValue *grid, size_t level, double epsilon, bool biharmonic) {
  size_t n = (size_t)pow(2, level);
  if (level > MIN_LEVEL) {
    /* Generate a coarser grid and solve that first to get good starting values */
    size_t level1 = level - 1, n1 = (size_t)pow(2, level1);
    struct GridValue *grid1 = (struct GridValue *)malloc(n1 * n1 * sizeof(struct GridValue));
    for (size_t i = 0; i < n1; ++i)
      for (size_t j = 0; j < n1; ++j) {
        grid1[j*n1+i].value = 0.0;
        int count = 0;
        if (grid[2*j*n+2*i].boundary) {
          ++count;
          grid1[j*n1+i].value += grid[2*j*n+2*i].value;
        }
        if (grid[2*j*n+2*i+1].boundary) {
          ++count;
          grid1[j*n1+i].value += grid[2*j*n+2*i+1].value;
        }
        if (grid[(2*j+1)*n+2*i].boundary) {
          ++count;
          grid1[j*n1+i].value += grid[(2*j+1)*n+2*i].value;
        }
        if (grid[(2*j+1)*n+2*i+1].boundary) {
          ++count;
          grid1[j*n1+i].value += grid[(2*j+1)*n+2*i+1].value;
        }
        if (count > 0) {
          grid1[j*n1+i].boundary = true;
          grid1[j*n1+i].value /= (double)count;
        }
        else
          grid1[j*n1+i].boundary = false;
      }
    solve(grid1, level1, epsilon, biharmonic);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        if (!grid[j*n+i].boundary)
          grid[j*n+i].value = grid1[(j/2)*n1+i/2].value;
    free(grid1);
  }

  /* Solve by iteration */
  if (biharmonic)
    solveBiharmonic(grid, n, epsilon);
  else
    solveHarmonic(grid, n, epsilon);
}

struct HarmonicMap *harmonic_create(double *min, double *max, size_t levels) {
  /* Create the grid */
  size_t n = (size_t)pow(2, levels);
  struct HarmonicMap *map = (struct HarmonicMap *)malloc(sizeof(struct HarmonicMap));
  map->levels = levels;
  map->size = n;
  map->grid = (struct GridValue *)malloc(n * n * sizeof(struct GridValue));

  /* Add a margin of 2.5% on all sides */
  double length = MAX(max[0] - min[0], max[1] - min[1]);
  map->offset[0] = min[0] - length * 0.025;
  map->offset[1] = min[1] - length * 0.025;
  map->scaling = (double)n / length / 1.05;

  /* Initialize cells */
  for (size_t i = 0; i < n * n; ++i) {
    map->grid[i].boundary = false;
    map->grid[i].value = 0.0;
  }

  return map;
}

void harmonic_add_point(struct HarmonicMap *map, double *point) {
  size_t n = map->size;
  int x = (int)round((point[0] - map->offset[0]) * map->scaling);
  int y = (int)round((point[1] - map->offset[1]) * map->scaling);
  map->grid[y*n+x].boundary = true;
  map->grid[y*n+x].value = point[2];
}

void harmonic_add_line(struct HarmonicMap *map, double *from, double *to) {
  size_t n = map->size;
  int x0 = (int)round((from[0] - map->offset[0]) * map->scaling);
  int y0 = (int)round((from[1] - map->offset[1]) * map->scaling);
  double v0 = from[2];
  int x1 = (int)round((to[0] - map->offset[0]) * map->scaling);
  int y1 = (int)round((to[1] - map->offset[1]) * map->scaling);
  double v1 = to[2];
  int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = (dx > dy ? dx : -dy) / 2, e2;
  while (true) {
    double ratio;             /* linear interpolation along the sides */
    if (err > 0)
      ratio = (double)abs(x1 - x0) / (double)dx;
    else
      ratio = (double)abs(y1 - y0) / (double)dy;
    map->grid[y0*n+x0].boundary = true;
    map->grid[y0*n+x0].value = v0 * ratio + v1 * (1.0 - ratio);
    if (x0 == x1 && y0 == y1) break;
    e2 = err;
    if (e2 > -dx) { err -= dy; x0 += sx; }
    if (e2 <  dy) { err += dx; y0 += sy; }
  }
}

void harmonic_solve(struct HarmonicMap *map, double epsilon, bool biharmonic) {
  solve(map->grid, map->levels, epsilon, biharmonic);
}

bool inside_map(struct HarmonicMap *map, size_t i, size_t j) {
  return i >= 0 && j >= 0 && i < map->size && j < map->size;
}

bool harmonic_eval(struct HarmonicMap *map, double *point, double *value) {
  size_t n = map->size;
  double x = (point[0] - map->offset[0]) * map->scaling;
  double y = (point[1] - map->offset[1]) * map->scaling;
  int i = (int)round(x), j = (int)round(y);
  
  if (!(inside_map(map, i, j) ||
        inside_map(map, i, j + 1) ||
        inside_map(map, i + 1, j) ||
        inside_map(map, i + 1, j + 1)))
    return false;               /* The point is outside the region */

  *value = map->grid[j*n+i].value * (1.0 - y + j) * (1.0 - x + i);
  *value += map->grid[(j+1)*n+i].value * (y - j) * (1.0 - x + i);
  *value += map->grid[j*n+i+1].value * (1.0 - y + j) * (x - i);
  *value += map->grid[(j+1)*n+i+1].value * (y - j) * (x - i);

  return true;
}

void harmonic_write_ppm(struct HarmonicMap *map, char *filename) {
  size_t n = map->size;
  FILE *f = fopen(filename, "w");
  fprintf(f, "P3\n%zu %zu\n255\n", n, n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j)
      if (map->grid[j*n+i].boundary)
        fprintf(f, "255 0 0 ");
      else
        fprintf(f, "0 0 %d ", (int)round(map->grid[j*n+i].value * 255.0));
    fprintf(f, "\n");
  }
  fclose(f);
}

void harmonic_free(struct HarmonicMap *map) {
  free(map->grid);
  free(map);
}
