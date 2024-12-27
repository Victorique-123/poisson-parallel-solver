# 泊松方程并行求解

莫斯科国立大学计算数学与控制论系  
超级计算机建模和技术课程

[中文](./README.md) | [Русский](./README.ru.md)


### 项目简介
本项目实现了使用不同并行计算方法求解泊松方程的程序，包括顺序执行、OpenMP、MPI以及MPI+OpenMP混合方式。采用有限差分法进行离散化，并使用虚拟区域法处理曲线边界。

### 文件说明
* `sequential.cpp` - 基础版本，无并行化，使用共轭梯度法求解
* `openmp.cpp` - OpenMP版本，使用最速下降法
* `openmp_cg.cpp` - OpenMP版本，使用共轭梯度法
* `mpi.cpp` - MPI版本程序
* `mpi_omp.cpp` - MPI+OpenMP混合版本程序
* `run_tests.sh` - 用于自动化测试的bash脚本
* `log/` - 包含实验过程中生成的日志文件

### 编译方法

```bash
# 顺序版本
g++ -fopenmp sequential.cpp -O3 -std=c++11 -o sequential

# OpenMP版本
g++ -fopenmp openmp.cpp -O3 -std=c++11 -o openmp

# MPI版本（标准环境）
mpicxx mpi.cpp -O3 -std=c++11 -o mpi

# MPI版本（在Polus上）
module load SpectrumMPI/10.1.0
mpixlC -std=c++11 -O3 -qstrict mpi.cpp -o mpi

# MPI+OpenMP版本（标准环境）
mpicxx -fopenmp mpi_omp.cpp -O3 -std=c++11 -o mpi_omp

# MPI+OpenMP版本（在Polus上）
module load SpectrumMPI/10.1.0
module load OpenMPI/4.0.2
mpicxx -fopenmp mpi_omp.cpp -O3 -std=c++11 -o mpi_omp
```

### 运行方法

```bash
# 顺序版本
./sequential {M} {N} {求解精度}

# OpenMP版本
./openmp {M} {N} {求解精度} {OpenMP线程数}

# MPI版本
mpirun -np {MPI进程数} ./mpi {M} {N} {求解精度}

# MPI+OpenMP版本
mpirun -np {MPI进程数} ./mpi_omp {M} {N} {求解精度} {每个进程的OpenMP线程数}
```
