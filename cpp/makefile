CXX = g++
CXXFLAGS = -std=c++11 -I/usr/lib/x86_64-linux-gnu/ -I/home/ubuntu/workspace/DateTime/Code/
INCLUDES = 
LD = icc
LDFLAGS := -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_program_options -lboost_system -L/home/ubuntu/workspace/DateTime/Code/

#/home/ubuntu/workspace/boost
#/home/ubuntu/workspace/boost/stage/lib

SOURCES = matrix.cpp datetime.cpp
HEADERS = matrix.h datetime.h
OBJ     = main.o matrix.o datetime.o

.SUFFIXES: .cpp

all: mat

mat: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o mat $(LDFLAGS)
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
	
matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp

datetime.o: ~/workspace/DateTime/Code/datetime.cpp ~/workspace/DateTime/Code/datetime.h
	$(CXX) $(CXXFLAGS) -c ~/workspace/DateTime/Code/datetime.cpp

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.lst
	rm -f output.*