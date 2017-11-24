CC = g++
OBJECTS = main.o graph.o asa103.o gpb.o utils.o
LDFLAGS = -larmadillo
CPPFLAGS = -std=c++11 -O2

gpb : $(OBJECTS)
	$(CC) -o gpb $(OBJECTS) $(LDFLAGS) $(CPPFLAGS)


clean:
	rm -rf *.o
