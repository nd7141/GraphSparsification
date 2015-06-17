all: Harvester

clean:
	-rm Harvester.o Harvester Harvester.d

Harvester.o: main.cpp
	g++ -O3 -g3 -Wall -std=c++11 -c -fmessage-length=0 -MMD -MP -MF"Harvester.d" -MT"Harvester.d" -o "Harvester.o" "./main.cpp"

Harvester: Harvester.o
	g++  -o "Harvester"  ./Harvester.o 	