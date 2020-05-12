all:
	icc -Ofast -Wall -std=c++11 mimd_interpreter-papi.cpp -o mimd_interpreter-transp-papi.x -funroll-loops -march=native -mtune=native -xavx2  -m64 -fno-alias -L/home/fazenda/papi-install/lib -lpapi -I/home/fazenda/papi-install/include
	icc -Ofast -Wall -std=c++11 interpreter_transpose-papi.cpp -o interpreter_transpose-transp-papi.x -funroll-loops -march=native -mtune=native -xavx2  -m64 -fno-alias -L/home/fazenda/papi-install/lib -lpapi -I/home/fazenda/papi-install/include
