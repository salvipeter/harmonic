all: libharmonic.a test

CFLAGS=-std=c99 -Wall -pedantic
TARGET_DIR=../target

libharmonic.a: $(TARGET_DIR)/harmonic.o
	ar rcs $(TARGET_DIR)/$@ $<

harmonic.o: harmonic.c harmonic.h
	gcc $(CFLAGS) -c -o $(TARGET_DIR)/$@ $<

test: test.c $(TARGET_DIR)/libharmonic.a
	gcc $(CFLAGS) -o $(TARGET_DIR)/$@ $< -L$(TARGET_DIR) -lharmonic -lm
