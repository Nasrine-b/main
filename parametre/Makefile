CC=gcc
CCFLAGS= -Wall -Werror -std=c99 -g 
LIBFLAGS= -lm
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
EXEC= exe

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $^ -o $@ $(LIBFLAGS)  

%.o: %.c
	$(CC) $(CCFLAGS) -o $@ -c $<

.depends:
	gcc -MM $(SRC) > .depends
-include .depends

clean:
	rm -f $(OBJ) $(EXEC)

distclean : clean