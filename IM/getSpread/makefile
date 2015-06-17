all: runCascade

clean:
	-rm runCascade.o runCascade runCascade.d

runCascade.o: runCascade.cpp
	g++ -O3 -g3 -Wall -std=c++11 -c -fmessage-length=0 -MMD -MP -MF"runCascade.d" -MT"runCascade.d" -o "runCascade.o" "./runCascade.cpp"

runCascade: runCascade.o
	g++  -o "runCascade"  ./runCascade.o 	