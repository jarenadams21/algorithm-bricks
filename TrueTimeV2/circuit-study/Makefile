CC = clang
CFLAGS = -std=c11 -O2 -Wall
LDLIBS = -lm

# Automatically gather all .c files in this directory
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# Default target
all: quantum_sdk_c
.PHONY: all clean

quantum_sdk_c: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDLIBS)

# Explicit rule for main.o in case pattern matching fails
main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

# Generic rule for any .c to .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) quantum_sdk_c