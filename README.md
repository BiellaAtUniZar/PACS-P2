# Lab 2

## Build project

``` 
mkdir build
cd build/
cmake ..
make -j4
```

## Run scripts

Run standard multiplication:
```
cd build/
./standar_multi arg1 arg2
```

where arg1 is the size N of the matrix and arg2 is the number of iterations to compute the mean execution time.


Run eigen multiplication:
```
cd build/
./eigen_multi arg1 arg2
```

where arg1 is the size N of the matrix and arg2 is the number of iterations to compute the mean execution time.
