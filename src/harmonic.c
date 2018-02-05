#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

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

#define MIN_LEVEL 3

#define MAX(a,b) (((a)>(b))?(a):(b))
#define TRUE 1
#define FALSE 0
typedef char Bool;

typedef struct PointStack_ {
  unsigned int x, y;
  struct PointStack_ *next;
} PointStack;

PointStack *push(PointStack *ps, unsigned int x, unsigned int y) {
  PointStack *e = (PointStack *)malloc(sizeof(PointStack));
  e->x = x;
  e->y = y;
  e->next = ps;
  return e;
}

PointStack *pop(PointStack *ps) {
  PointStack *next = ps->next;
  free(ps);
  return next;
}

void flood_fill(struct HarmonicMap *map, unsigned int x, unsigned int y) {
  unsigned int n = map->size;
  PointStack *ps = push(NULL, x, y);
  do {
    x = ps->x; y = ps->y;
    ps = pop(ps);
    if (map->grid[y*n+x].type != UNTYPED)
      continue;
    map->grid[y*n+x].type = EXTERIOR;
    if (x > 0)
      ps = push(ps, x - 1, y);
    if (x < n - 1)
      ps = push(ps, x + 1, y);
    if (y > 0)
      ps = push(ps, x, y - 1);
    if (y < n - 1)
      ps = push(ps, x, y + 1);
  } while (ps);
}

void solve(struct GridValue *grid, unsigned int level, double epsilon) {
  unsigned int n = (unsigned int)pow(2, level);
  if (level > MIN_LEVEL) {
    /* Generate a coarser grid and solve that first to get good starting values */
    unsigned int level1 = level - 1, n1 = (unsigned int)pow(2, level1);
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
  do {
    change = 0.0;
    unsigned int count = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j)
        if (grid[j*n+i].type == INTERIOR) {
          int neighbors = 0;
          double old_value = grid[j*n+i].value;
          grid[j*n+i].value = 0.0;
          if (j > 0 && grid[(j-1)*n+i].type != EXTERIOR) {
            ++neighbors;
            grid[j*n+i].value += grid[(j-1)*n+i].value;
          }
          if (i > 0 && grid[j*n+i-1].type != EXTERIOR) {
            ++neighbors;
            grid[j*n+i].value += grid[j*n+i-1].value;
          }
          if (j < n && grid[(j+1)*n+i].type != EXTERIOR) {
            ++neighbors;
            grid[j*n+i].value += grid[(j+1)*n+i].value;
          }
          if (i < n && grid[j*n+i+1].type != EXTERIOR) {
            ++neighbors;
            grid[j*n+i].value += grid[j*n+i+1].value;
          }
          assert(neighbors > 0);
          grid[j*n+i].value /= (double)neighbors;
          ++count;
          change += fabs(grid[j*n+i].value - old_value);
        }
    if (count)
      change /= (double)count;
  } while (change > epsilon);
}

struct HarmonicMap *harmonic_init(unsigned int size, double *points, unsigned int levels,
                                  double epsilon) {
  /* Create the grid */
  unsigned int n = (unsigned int)pow(2, levels);
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
  map->scaling = (double)n / length / 1.05;

  /* Fill boundary cells */
  for (unsigned int i = 0; i < n * n; ++i)
    map->grid[i].type = UNTYPED;
  int x0 = (int)round((points[size*3-3] - map->offset[0]) * map->scaling);
  int y0 = (int)round((points[size*3-2] - map->offset[1]) * map->scaling);
  double v0 = points[size*3-1];
  for (unsigned int i = 0; i < size; ++i) {
    /* line drawing from Rosetta Code [Bresenham's algorithm] */
    int x1 = (int)round((points[i*3] - map->offset[0]) * map->scaling);
    int y1 = (int)round((points[i*3+1] - map->offset[1]) * map->scaling);
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
  flood_fill(map, 0, 0);
  for (unsigned int i = 0; i < n * n; ++i)
    if (map->grid[i].type == UNTYPED) {
      map->grid[i].type = INTERIOR;
      map->grid[i].value = 0.0;
    }

  /* Add bogus values to exterior cells near the boundary to help evaluation */
  for (unsigned int y = 1; y < n - 1; ++y)
    for (unsigned int x = 1; x < n - 1; ++x) {
      unsigned int i = y * n + x;
      if (map->grid[i].type == EXTERIOR) {
        unsigned int count = 0;
        double value = 0.0;
        for (int dy = -1; dy <= 1; ++dy) {
          unsigned int j = i + dy * n;
          for (int dx = -1; dx <= 1; ++dx) {
            if (map->grid[j+dx].type == BOUNDARY) {
              ++count;
              value += map->grid[j+dx].value;
            }
          }
        }
        if (count > 0)
          map->grid[i].value = value / count;
      }
    }

  /* Find the solution for the discrete problem */
  solve(map->grid, levels, epsilon);

  return map;
}

bool inside_map(struct HarmonicMap *map, unsigned int i, unsigned int j) {
  return i >= 0 && j >= 0 && i < map->size && j < map->size;
}

bool harmonic_eval(struct HarmonicMap *map, double *point, double *value) {
  unsigned int n = map->size;
  double x = (point[0] - map->offset[0]) * map->scaling;
  double y = (point[1] - map->offset[1]) * map->scaling;
  int i = (int)round(x), j = (int)round(y);
  
  if (!((inside_map(map, i, j) && map->grid[j*n+i].type != EXTERIOR) ||
        (inside_map(map, i, j + 1) && map->grid[(j+1)*n+i].type != EXTERIOR) ||
        (inside_map(map, i + 1, j) && map->grid[j*n+i+1].type != EXTERIOR) ||
        (inside_map(map, i + 1, j + 1) && map->grid[(j+1)*n+i+1].type != EXTERIOR)))
    return false;               /* The point is outside the region */

  *value = map->grid[j*n+i].value * (1.0 - y + j) * (1.0 - x + i);
  *value += map->grid[(j+1)*n+i].value * (y - j) * (1.0 - x + i);
  *value += map->grid[j*n+i+1].value * (1.0 - y + j) * (x - i);
  *value += map->grid[(j+1)*n+i+1].value * (y - j) * (x - i);

  return true;
}

void harmonic_write_ppm(struct HarmonicMap *map, char *filename) {
  unsigned int n = map->size;
  FILE *f = fopen(filename, "w");
  fprintf(f, "P3\n%d %d\n255\n", n, n);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j)
      if (map->grid[j*n+i].type == EXTERIOR)
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

unsigned int harmonic_mesh_size(struct HarmonicMap *map, unsigned int downsampling) {
  unsigned int n = map->size, count = 0;
  unsigned int k = (unsigned int)pow(2, downsampling);
  for (unsigned int j = k; j < n; j += k)
    for (unsigned int i = 0; i < n - k; i += k)
      if (map->grid[j*n+i].type != EXTERIOR)
        ++count;
  return count;
}

unsigned int harmonic_mesh(struct HarmonicMap *map, unsigned int downsampling,
                           double *vertices, unsigned int *triangles) {
  unsigned int n = map->size, index = 0, count = 0;
  unsigned int k = (unsigned int)pow(2, downsampling);
  unsigned int *row = (unsigned int *)malloc(sizeof(unsigned int) * n);

  /* Note: because of the enlargement, there can be no vertices in the first row,
     or at the end of a row, thus simplifying the algorithm. */

  for (unsigned int j = k; j < n; j += k) {
    for (unsigned int i = 0; i < n - k; i += k) {
      if (map->grid[j*n+i].type == EXTERIOR)
        continue;

      if (map->grid[(j-k)*n+i+k].type == EXTERIOR) {
        /* no NE */
        if (map->grid[(j-k)*n+i].type != EXTERIOR &&
            map->grid[j*n+i+k].type != EXTERIOR) {
          /* N & E */
          triangles[3*count+0] = index;
          triangles[3*count+1] = row[i];
          triangles[3*count+2] = index + 1;
          ++count;
        }
      } else {
        if (map->grid[(j-k)*n+i].type != EXTERIOR) {
          /* N & NE */
          triangles[3*count+0] = index;
          triangles[3*count+1] = row[i];
          triangles[3*count+2] = row[i+k];
          ++count;
        }
        if (map->grid[j*n+i+k].type != EXTERIOR) {
          /* E & NE */
          triangles[3*count+0] = index;
          triangles[3*count+1] = row[i+k];
          triangles[3*count+2] = index + 1;
          ++count;
        }
      }

      /* Update row */
      vertices[2*index+0] = (double)i / map->scaling + map->offset[0];
      vertices[2*index+1] = (double)j / map->scaling + map->offset[1];
      row[i] = index++;
    }
  }

  free(row);

  return count;
}
