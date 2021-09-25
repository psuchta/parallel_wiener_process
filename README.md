# parallel_wiener_process

## How to compile OpenMP program
*macOS*
```
g++ -Xpreprocessor -fopenmp -std=c++11 wiener.cpp -o wiener -lomp
```

## How to compile Open MPI program
*macOS*
```
mpic++ -std=c++11 -o wiener ./wiener.cpp

mpirun -np 4 ./wiener