#!/bin/bash

# 加载必要的模块
module load SpectrumMPI/10.1.0
module load OpenMPI/4.0.2

# 编译MPI程序
compile_programs() {
    echo "Compiling MPI+OpenMP program..."
    mpicxx -fopenmp -O3 -o mpi_opm mpi_opm.cpp -std=c++11
    
    if [ $? -eq 0 ]; then
        echo "Compilation successful"
    else
        echo "Compilation failed"
        exit 1
    fi
}

# 提交MPI+OpenMP作业 (80x90) 1
submit_mpi_opm_jobs_80x90_first() {
    echo "Submitting MPI jobs for 80x90 grid..."
        bsub -n 1 -q normal \
             -o mpi_opm_80x90.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_opm_80x90_${procs}p \
             "ulimit -s 102400; mpirun -np 1 ./mpi_opm 80 90 1.00e-08 1"
}

# 提交MPI+OpenMP作业 (80x90) 1
submit_mpi_opm_jobs_160x180_first() {
    echo "Submitting MPI jobs for 80x90 grid..."
        bsub -n 1 -q normal \
             -o mpi_opm_160x180.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_opm_160x180_${procs}p \
             "ulimit -s 102400; mpirun -np 1 ./mpi_opm 160 180 2.50e-09 1"
}


# 提交MPI+OpenMP作业 (80x90)
submit_mpi_opm_jobs_80x90() {
    echo "Submitting MPI jobs for 80x90 grid..."
    
    for procs in 1 2 4 8; do
        bsub -n 2 -q normal \
             -o mpi_opm_80x90.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_opm_80x90_${procs}p \
             "ulimit -s 102400; mpirun -np 2 ./mpi_opm 80 90 1.00e-08 ${procs}"
    done
}

# 提交MPI+OpenMP作业 (160x180)
submit_mpi_opm_jobs_160x180() {
    echo "Submitting MPI jobs for 160x180 grid..."
    
    for procs in 1 2 4 8; do
        bsub -n 4 -q normal \
             -o mpi_opm_160x180.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_opm_160x180_${procs}p \
             "ulimit -s 102400; mpirun -np 4 ./mpi_opm 160 180 2.50e-09 ${procs}"
    done
}

# 主函数
main() {
    compile_programs
    
    echo "Do you want to submit all jobs? [y/n]"
    read answer
    
    if [ "$answer" = "y" ]; then
        submit_mpi_opm_jobs_80x90_first
        submit_mpi_opm_jobs_160x180_first
        submit_mpi_opm_jobs_80x90
        submit_mpi_opm_jobs_160x180
        echo "All jobs submitted successfully"
    else
        echo "Available options:"
        echo "1: Submit MPI+OpenMP jobs (80x90_first)"
        echo "2: Submit MPI+OpenMP jobs (160x180_first)"
        echo "3: Submit MPI+OpenMP jobs (80x90)"
        echo "4: Submit MPI+OpenMP jobs (160x180)"
        echo "Enter option number (or 'q' to quit):"
        
        read option
        case $option in
            1) submit_mpi_opm_jobs_80x90_first ;;
            2) submit_mpi_opm_jobs_160x180_first ;;
            3) submit_mpi_opm_jobs_80x90 ;;
            4) submit_mpi_opm_jobs_160x180 ;;
            q) echo "Exiting..." ;;
            *) echo "Invalid option" ;;
        esac
    fi
}

# 执行主函数
main