# ACSE 6  Solving the Wave Equation

[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

The project is a C++ package to implement to solve the Wave Equation in parallel. The python file is used to collect the .dat files and output pictures and video.

## Table of Contents

- [About project](#about-project)
- [Configuration and Installation](#configuration-and-installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Reference](#reference)
- [License](#license)


## About project

The aim of this assignment is to write a parallel solver for the wave equation. The boundary conditions include Dirichlet boundary and Neumann boundary. The decomposition included stripe decomposition and grid decomposition.

After getting the solution, the post-processing module can convert the .dat file into images and videos.


## Configuration and installation

To run this project, you should download some files, which list below

1. MPI : download the MPI and configure it in the visual studio. The URL is

```
http://www.mpich.org/downloads/
```

2. download the numpy and Image module

```
pip install numpy
pip install Image
```


## Usage

This section will introduce on how to run this project and give a few simple examples.

### Setting
The default size of the array is 100*100. You can change the size through changing imax and jmax in `Matrix.h`. In addition, some other parameters like dt_out, x_max, y_max also can be changed.


The boundary conditions can be changed in `main.cpp`. 

In addition, before running the code, you must create four folders whose name are  `out`, `Pictures`, `result` and `test`. Please don't change the files provided in `test` folder because these files are used in testing. The other three folder is used to  receive .dat files or .png files.

### Running
Run the main file in cmd or powershell with this command:

```
mpiexec -n 10 ACSE6.2.exe
```

10 is the number of the processes. You can change other integers than more than 0.  ACSE6.2 is my project name. When you change the Project name, yous should input the name coincide your project name

If you don't change the code and run the code directly, you will get

```
p0: 0 9
p1: 9 9
p2: 18 10
p3: 28 10
p4: 38 10
p5: 48 10
p6: 58 10
p7: 68 10
p8: 78 10
p9: 88 10
output: 300     t: 30   iteration: 2969
Divide 10 into 2 by 5 grid
p: 0  row_start 0  col start 0 row_chunk: 49 col_chunk: 19
p: 1  row_start 0  col start 19 row_chunk: 49 col_chunk: 19
p: 2  row_start 0  col start 38 row_chunk: 49 col_chunk: 20
p: 3  row_start 0  col start 58 row_chunk: 49 col_chunk: 20
p: 4  row_start 0  col start 78 row_chunk: 49 col_chunk: 20
p: 5  row_start 49  col start 0 row_chunk: 49 col_chunk: 19
p: 6  row_start 49  col start 19 row_chunk: 49 col_chunk: 19
p: 7  row_start 49  col start 38 row_chunk: 49 col_chunk: 20
p: 8  row_start 49  col start 58 row_chunk: 49 col_chunk: 20
p: 9  row_start 49  col start 78 row_chunk: 49 col_chunk: 20
output: 300     t: 30   iteration: 2969
```

This sizes of the processes is 10 in this example and each solver generate 300 .dat files in `out` folder .

### Collect file
Open the Collect.py. You will find some parameters. The default parameters are


```
dir_out = os.getcwd() + '\out'
dir_test = os.getcwd() + '\\test'
dir_result = os.getcwd() + '\\result'

p = 10
cnt = 300
imax = 100
jmax = 100

fps = 24
decomposition = 0
```

p is the size of the processes that you input in cmd or powershell. The result output = 300, therefore, set cnt = 300. imax and jmax are same to the imax and jamx in `Matrix.h`. fps is the default value.

After running the function `collectgrid` or `collectstrip`, the .dat file will be  generated in `result` folder. 

### Visualisation

Visualisation relies on the files containing the value information in `result` folder. 

After reading in this file, the function `outPictures` will generate .png files in `Pictures` folder. Then, using the function `imageTovideo` to collect the pictures to generate the vedio.



## Documentation

To get documentation of this project, open `index.html` in `documentation/html` after installation.


## Reference
ACSE 6 Lecture 2 Exercise 3

ACSE 6 Lecture 4 Exercise 2

ACSE 6 Lecture 6 Exercise 1

## License

[MIT](LICENSE) © acse-2020
