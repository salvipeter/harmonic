all: libharmonic.so test

CFLAGS=-std=c99 -Wall -pedantic
TARGET_DIR=../target

libharmonic.so: harmonic.c point_stack.c
	gcc -O3 -fpic -shared $(CFLAGS) -o $(TARGET_DIR)/$@ $^

test: test.c libharmonic.so
	gcc $(CFLAGS) -o $(TARGET_DIR)/$@ $< -L$(TARGET_DIR) -lharmonic -lm
