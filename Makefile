CC = g++
SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)
LDFLAGS = -larmadillo
CPPFLAGS = -std=c++11 -O2

gpb : $(OBJS)
	$(CC) -o gpb $(OBJS) $(LDFLAGS) $(CPPFLAGS)

clean:
	rm -rf src/*.o
