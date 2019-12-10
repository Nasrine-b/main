CC=gcc
CCFLAGS= -Wall -g 
LIBFLAGS= -lm
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
EXEC= exe


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LIBFLAGS) $^ -o $@  

%.o: %.cc
	$(CC) $(CCFLAGS) -o $@ -c $<

.depend:
	g++ -MM $(SRC) > .depends
-include .depends

clean:
	rm -f $(OBJ) $(EXEC)

distclean : clean

