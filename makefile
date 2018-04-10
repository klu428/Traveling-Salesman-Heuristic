CXX = g++
RM=rm -f
CXXFLAGS = -std=c++0x
SRCS = tsp.cpp
HEADERS =
OBJS = *.o
OUT = *.tour
UNIX = circular
WINDOWS = *.exe

tsp: ${SRCS} ${HEADERS}
	${CXX} ${CXXFLAGS} ${SRCS} -o tsp

#put the directory back to its starting state (just .cpp, .hpp, and makefile files)
clean:
	${RM} ${OBJS} ${OUT} ${UNIX} ${WINDOWS}
