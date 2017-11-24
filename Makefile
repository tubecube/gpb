CC = g++
SRCS = $(wildcard src/utils.cpp src/graph.cpp src/gpb.cpp src/main.cpp src/asa103.cpp)
OBJS = $(SRCS:.cpp = .o)
LDFLAGS = -larmadillo
CPPFLAGS = -std=c++11 -O2

gpb : $(OBJS)
	$(CC) -o gpb $(OBJS) $(LDFLAGS) $(CPPFLAGS)

clean:
	rm -rf *.o
