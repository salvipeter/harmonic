#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

#define MIN_LEVEL 3

#define MAX(a,b) (((a)>(b))?(a):(b))
#define TRUE 1
#define FALSE 0
typedef char Bool;

#define RECURSIVE_FLOOD_FILL
#ifdef RECURSIVE_FLOOD_FILL

void flood_fill(struct HarmonicMap *map, unsigned int x, unsigned int y) {
  unsigned int n = map->size;
  if (map->grid[y*n+x].type != UNTYPED)
    return;
  map->grid[y*n+x].type = EXTERIOR;
  if (x > 0)
    flood_fill(map, x - 1, y);
  if (x < n - 1)
    flood_fill(map, x + 1, y);
  if (y > 0)
    flood_fill(map, x, y - 1);
  if (y < n - 1)
    flood_fill(map, x, y + 1);
}

#else // !RECURSIVE_FLOOD_FILL

/* Flood fill implementation based on the fixed-memory algorithm pseudocode on Wikipedia */

struct Position {
  Bool valid;
  int x, y;
  enum Direction { LEFT, RIGHT, UP, DOWN } dir;
};

Bool is_valid(struct HarmonicMap *map, struct Position *p) {
  unsigned int n = map->size;
  if (p->x < 0 || p->y < 0 || p->x >= n || p->y >= n)
    return FALSE;
  return TRUE;
}

struct GridValue *get_cell(struct HarmonicMap *map, struct Position *p) {
  return &map->grid[p->y*map->size+p->x];
}

void move_forward(struct Position *p) {
  switch (p->dir) {
  case LEFT:  p->x--; break;
  case RIGHT: p->x++; break;
  case UP:    p->y--; break;
  case DOWN:  p->y++;
  }
}

void turn_left(struct Position *p) {
  switch (p->dir) {
  case LEFT:  p->dir = DOWN;  break;
  case RIGHT: p->dir = UP;    break;
  case UP:    p->dir = LEFT;  break;
  case DOWN:  p->dir = RIGHT; break;
  }
}

void turn_right(struct Position *p) {
  switch (p->dir) {
  case LEFT:  p->dir = UP;    break;
  case RIGHT: p->dir = DOWN;  break;
  case UP:    p->dir = RIGHT; break;
  case DOWN:  p->dir = LEFT;  break;
  }
}

void turn_around(struct Position *p) {
  turn_right(p);
  turn_right(p);
}

Bool front_cell_empty(struct HarmonicMap *map, struct Position *p) {
  struct Position q = *p;
  move_forward(&q);
  return is_valid(map, &q) && get_cell(map, &q)->type == UNTYPED;
}

Bool right_cell_empty(struct HarmonicMap *map, struct Position *p) {
  struct Position q = *p;
  turn_right(&q);
  move_forward(&q);
  return is_valid(map, &q) && get_cell(map, &q)->type == UNTYPED;
}

Bool left_cell_empty(struct HarmonicMap *map, struct Position *p) {
  struct Position q = *p;
  turn_left(&q);
  move_forward(&q);
  return is_valid(map, &q) && get_cell(map, &q)->type == UNTYPED;
}

Bool back_cell_empty(struct HarmonicMap *map, struct Position *p) {
  struct Position q = *p;
  turn_around(&q);
  move_forward(&q);
  return is_valid(map, &q) && get_cell(map, &q)->type == UNTYPED;
}

Bool front_left_cell_empty(struct HarmonicMap *map, struct Position *p) {
  struct Position q = *p;
  move_forward(&q);
  turn_left(&q);
  move_forward(&q);
  return is_valid(map, &q) && get_cell(map, &q)->type == UNTYPED;
}

Bool back_left_cell_empty(struct HarmonicMap *map, struct Position *p) {
  struct Position q = *p;
  turn_left(&q);
  move_forward(&q);
  turn_left(&q);
  move_forward(&q);
  return is_valid(map, &q) && get_cell(map, &q)->type == UNTYPED;
}

void flood_fill(struct HarmonicMap *map, unsigned int x, unsigned int y) {
  /* this algorithm goes along the right-hand walls like in a maze
     mark is used to save a two-way position that maybe we have no path back to
     - the mark can be cleared if we find a way back to it
     - if we visit it next time in the same direction, it can be painted
     - if in the other direction, there is a loop further down the path
       - backtrack shows that we are looking for this loop
       - mark2 is used to save a two-way position inside this loop
       - findloop shows if there was a fork in the road since backtracking started */
  struct Position cur, mark, mark2;
  Bool backtrack, findloop;
  unsigned char count;
  cur.valid = TRUE; cur.x = x; cur.y = y; cur.dir = LEFT;
  mark.valid = FALSE; mark2.valid = FALSE;
  backtrack = FALSE; findloop = FALSE;

  /* Find a wall */
  while (front_cell_empty(map, &cur)) {
    move_forward(&cur);
  }

  goto start;

 main_loop:
  /* Move forward and turn right if there is a path;
     set findloop if there is an alternative path and backtrack is set */
  move_forward(&cur);
  if (right_cell_empty(map, &cur)) {
    if (backtrack && !findloop && (front_cell_empty(map, &cur) || left_cell_empty(map, &cur)))
      findloop = TRUE;
    turn_right(&cur);
 paint:
    move_forward(&cur);
  }

 start:
  /* Count surrounding walls */
  count = 0;
  if (!front_cell_empty(map, &cur)) ++count;
  if (!right_cell_empty(map, &cur)) ++count;
  if (!left_cell_empty(map, &cur))  ++count;
  if (!back_cell_empty(map, &cur))  ++count;

  /* Turn until there is a wall on the right-hand side and an empty cell in front */
  if (count < 4) {
    do {
      turn_right(&cur);
    } while(front_cell_empty(map, &cur));
    do {
      turn_left(&cur);
    } while(!front_cell_empty(map, &cur));
  }

  switch (count) {
  case 1:
    /* If we are backtracking => set findloop and go along the wall
       If we are in findloop mode (and not backtracking) =>
         re-validate the mark, then go along the wall
       If both left corners are empty =>
         we can forget the mark, paint the cell and go forward
       If at least one of the left corners is painted =>
         we cannot paint, as this could subdivide the region -> go along the wall */
    if (backtrack)
      findloop = TRUE;
    else if (findloop) {
      if (!mark.valid)
        mark.valid = TRUE;
    } else if (front_left_cell_empty(map, &cur) && back_left_cell_empty(map, &cur)) {
      mark.valid = FALSE;
      get_cell(map, &cur)->type = EXTERIOR; /* Paint the current cell */
      goto paint;
    }
    break;
  case 2:
    if (!back_cell_empty(map, &cur)) {
      /* Exits: front & left
         If the front-left corner is painted =>
           we cannot paint, as this could subdivide the region -> go along the wall
         Otherwise we can forget the mark, paint the cell and go forward */
      if (front_left_cell_empty(map, &cur)) {
        mark.valid = FALSE;
        get_cell(map, &cur)->type = EXTERIOR; /* Paint the current cell */
        goto paint;
      }
    } else if (!mark.valid) {
      /* Exits: front & back, mark is not set =>
           we cannot paint, as this could cut us off -> save position and go along the wall */
      mark = cur;
      mark2.valid = FALSE;
      findloop = FALSE;
      backtrack = FALSE;
    } else {
      /* Exits: front & back, mark is set */
      if (!mark2.valid) {
        if (cur.x == mark.x && cur.y == mark.y) {
          /* We are at the (only) marked position */
          if (cur.dir == mark.dir) {
            /* Same direction => found a loop, so it is safe to paint
                 delete mark, and go the other way (exploring the left wall side) */
            mark.valid = FALSE;
            turn_around(&cur);
            get_cell(map, &cur)->type = EXTERIOR; /* Paint the current cell */
            goto paint;
          } else {
            /* There was a loop somewhere - go for another round with backtracking */
            backtrack = TRUE;
            findloop = FALSE;
            cur.dir = mark.dir;
          }
        } else if (findloop) {
          /* We are at an unmarked position wanting to find a loop =>
               place another mark and go along the wall */
          mark2 = cur;
        }
      } else {
        if (cur.x == mark.x && cur.y == mark.y) {
          /* We are back at mark, having explored the loop and marked a point of it =>
               go to mark2 - this is a point in the loop, so it is safe to paint
                 delete both marks, and go the other way (exploring the left wall side) */
          cur = mark2;
          mark.valid = FALSE;
          mark2.valid = FALSE;
          backtrack = FALSE;
          turn_around(&cur);
          get_cell(map, &cur)->type = EXTERIOR; /* Paint the current cell */
          goto paint;
        } else if (cur.x == mark2.x && cur.y == mark2.y) {
          /* We are back at mark2, so the loop is further down the path =>
               let this be the main mark, and continue in the old direction */
          mark = mark2;
          cur.dir = mark2.dir;
          mark2.valid = FALSE;
        }
      }
    }
    break;
  case 3:
    /* Forget mark, paint the current cell and go forward */
    mark.valid = FALSE;
    get_cell(map, &cur)->type = EXTERIOR; /* Paint the current cell */
    goto paint;
  case 4:
    /* Paint the last cell and exit */
    get_cell(map, &cur)->type = EXTERIOR;
    return;
  default:
    break;
  }

  goto main_loop;
}

#endif

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
    change /= (double)count;
  } while (change > epsilon);
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
  map->scaling = (double)n / length / 1.05;

  /* Fill boundary cells */
  for (unsigned int i = 0; i < n * n; ++i)
    map->grid[i].type = UNTYPED;
  int x0 = (points[size*3-3] - map->offset[0]) * map->scaling;
  int y0 = (points[size*3-2] - map->offset[1]) * map->scaling;
  double v0 = points[size*3-1];
  for (unsigned int i = 0; i < size; ++i) {
    /* line drawing from Rosetta Code [Bresenham's algorithm] */
    int x1 = (points[i*3] - map->offset[0]) * map->scaling;
    int y1 = (points[i*3+1] - map->offset[1]) * map->scaling;
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

  /* Find the solution for the discrete problem */
  solve(map->grid, levels, epsilon);

  return map;
}

double harmonic_eval(struct HarmonicMap *map, double *point) {
  unsigned int n = map->size;
  double x = (point[0] - map->offset[0]) * map->scaling;
  double y = (point[1] - map->offset[1]) * map->scaling;
  int i = (int)x, j = (int)y;
  if (i == map->size - 1)
    i = map->size - 2;
  if (j == map->size - 1)
    j = map->size - 2;
  /* TODO: check exterior cells */
  return map->grid[j*n+i].value * (1.0 - y + j) * (1.0 - x + i) +
    map->grid[(j+1)*n+i].value * (y - j) * (1.0 - x + i) +
    map->grid[j*n+i+1].value * (1.0 - y + j) * (x - i) +
    map->grid[(j+1)*n+i+1].value * (y - j) * (x - i);
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
        fprintf(f, "0 0 %d ", (int)(map->grid[j*n+i].value * 255.0));
    fprintf(f, "\n");
  }
  fclose(f);
}

void harmonic_free(struct HarmonicMap *map) {
  free(map->grid);
  free(map);
}
