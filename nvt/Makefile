
output: main.o particle.o energy.o box.o
	g++ main.o particle.o energy.o box.o -o output

main.o: main.cpp
	g++ -c main.cpp

particle.o: particle.cpp
	g++ -c particle.cpp

energy.o: energy.cpp
	g++ -c energy.cpp

box.o: box.cpp
	g++ -c box.cpp

clean: 
	rm *.o output

