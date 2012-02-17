# Imputation Makefile

CC = g++
CFLAGS = -g -W -Wall -Werror -Wno-unused
V = @

all: imputation

main.o: main.cpp
	@echo + cc main.cpp
	$(V)$(CC) $(CFLAGS) -c main.cpp

imputation.o: imputation.cpp
		@echo + cc imputation.cpp
		$(V)$(CC) $(CFLAGS) -c imputation.cpp

imputation: main.o imputation.o
	@echo + link main
	$(V)$(CC) $(CFLAGS) -o $@ main.o imputation.o

clean:
	@echo + clean
	$(V)rm -rf *.o imputation
