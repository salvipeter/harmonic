all: libharmonic.so test

CFLAGS=-std=c99 -Wall -pedantic -O3
TARGET_DIR=../target

libharmonic.so: harmonic.c
	gcc -fpic -shared $(CFLAGS) -o $(TARGET_DIR)/$@ $<

test: test.c libharmonic.so
	gcc $(CFLAGS) -o $(TARGET_DIR)/$@ $< -L$(TARGET_DIR) -lharmonic -lm
