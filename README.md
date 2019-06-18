# LinearSystems

This project contains methods for parsing a matrix using dense and CRS storage methods. Also, methods for solving a matrix using LU Factorization and Successive Over-Relaxation.

## Matrix Format
Matrix format must be exactly as the formats followed in the matrices from the folder `matrices`.


## Usage

### Compile the project
gcc *.c -o LinearSystems

### Arguments 
```
 --file [filename]                          The desired matrix 
 --lu                                       Method for LU factorization
 --sor  [omega] [tolerance] [iterations]    Method for successive over-relaxation
 --dense                                    For dense storage
 --crs                                      For compressed sparse rows storage
```

### Examples
#### Solving matrix using Dense Storage and LU method
```
./LinearSystems --file matrices/arc130.mtx --dense --lu 
```
#### Solving matrix using  Dense Storage and SOR method
![SOR method](https://wikimedia.org/api/rest_v1/media/math/render/svg/57760458e6cf8203f399cb09445a7989c0139cdc)

```
./LinearSystems --file matrices/arc130.mtx --dense --sor 0.5 0.00001 100
```
#### Solving matrix using CRS Storage and SOR method
```
./LinearSystems --file matrices/arc130.mtx --crs --sor 0.5 0.00001 100
```


