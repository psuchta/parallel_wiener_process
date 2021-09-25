# parallel_wiener_process

## Execution time

### Sequential
`10432.4 ms`

### Threads - 1 version
`7465.45 ms`
### Threads - 2 version
`7205.92 ms`
### OpenMP
`6782.73 ms`
### Open MPI
`2790.87 ms`


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