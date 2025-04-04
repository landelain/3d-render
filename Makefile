all : project1

project1 : project1.cpp
	g++ -fopenmp -o project1 project1.cpp

clean :
	rm -f project1