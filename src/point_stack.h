#pragma once

#include <stdbool.h>

struct point_stack_t *point_stack_new();

void point_stack_free(struct point_stack_t *ps);

void point_stack_push(struct point_stack_t *ps, size_t x, size_t y);

void point_stack_pop(struct point_stack_t *ps);

bool point_stack_empty(struct point_stack_t *ps);

void point_stack_top(struct point_stack_t *ps, size_t *x, size_t *y);
