#include <stdbool.h>
#include <stdlib.h>

typedef struct {
  size_t x, y;
} point_t;

typedef struct {
  size_t capacity, size;
  point_t *data;
} point_stack_t;

point_stack_t *point_stack_new() {
  point_stack_t *ps = (point_stack_t *)malloc(sizeof(point_stack_t));
  ps->capacity = 0;
  ps->size = 0;
  ps->data = NULL;
  return ps;
}

void point_stack_free(point_stack_t *ps) {
  free(ps->data);
  free(ps);
}

void point_stack_push(point_stack_t *ps, size_t x, size_t y) {
  if (ps->capacity == 0) {
    ps->capacity = 1;
    ps->data = (point_t *)malloc(sizeof(point_t));
  } else if (ps->size == ps->capacity) {
    ps->capacity *= 2;
    ps->data = (point_t *)realloc(ps->data, sizeof(point_t) * ps->capacity);
  }
  ps->data[ps->size].x = x;
  ps->data[ps->size].y = y;
  ps->size++;
}

void point_stack_pop(point_stack_t *ps) {
  ps->size--;
}

bool point_stack_empty(point_stack_t *ps) {
  return !ps->size;
}

void point_stack_top(point_stack_t *ps, size_t *x, size_t *y) {
  *x = ps->data[ps->size-1].x;
  *y = ps->data[ps->size-1].y;
}
