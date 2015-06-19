all: getPossibleWorlds

clean:
	-rm main.o getPossibleWorlds main.d

main.o: main.cpp
	g++ -O3 -g3 -Wall -std=c++11 -c -fmessage-length=0 -MMD -MP -MF"main.d" -MT"main.d" -o "main.o" "./main.cpp"

getPossibleWorlds: main.o
	g++  -o "getPossibleWorlds"  ./main.o 	