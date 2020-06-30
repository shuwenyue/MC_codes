
output: main.o particle.o energy.o box.o histogram.o
	g++ main.o particle.o energy.o box.o histogram.o -o output

main.o: main.cpp
	g++ -c main.cpp

particle.o: particle.cpp
	g++ -c particle.cpp

energy.o: energy.cpp
	g++ -c energy.cpp

box.o: box.cpp
	g++ -c box.cpp

histogram.o: histogram.cpp
	g++ -c histogram.cpp

clean: 
	rm *.o output

