#!/bin/bash

# 加载必要的模块
module load SpectrumMPI/10.1.0

# 编译MPI程序
compile_programs() {
    echo "Compiling MPI program..."
    mpixlC -std=c++11 -O3 -qstrict mpi.cpp -o mpi
    
    if [ $? -eq 0 ]; then
        echo "Compilation successful"
    else
        echo "Compilation failed"
        exit 1
    fi
}

# 提交MPI作业 (40x40)
submit_mpi_jobs() {
    echo "Submitting MPI jobs for 40x40 grid..."
    
    for procs in 1 2 4; do
        bsub -n ${procs} -q normal \
             -o mpi_40x40.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_40x40_${procs}p \
             "mpirun -np $procs ./mpi 40 40 1.00e-07"
    done
}

# 提交MPI作业 (80x90)
submit_mpi_jobs_80x90() {
    echo "Submitting MPI jobs for 80x90 grid..."
    
    for procs in 1 2 4; do
        bsub -n ${procs} -q normal \
             -o mpi_80x90.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_80x90_${procs}p \
             "mpirun -np $procs ./mpi 80 90 1.00e-08"
    done
}

# 提交MPI作业 (160x180)
submit_mpi_jobs_160x180() {
    echo "Submitting MPI jobs for 160x180 grid..."
    
    for procs in 1 2 4; do
        bsub -n ${procs} -q normal \
             -o mpi_160x180.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J mpi_160x180_${procs}p \
             "mpirun -np $procs ./mpi 160 180 2.50e-09"
    done
}

# 主函数
main() {
    compile_programs
    
    echo "Do you want to submit all jobs? [y/n]"
    read answer
    
    if [ "$answer" = "y" ]; then
        #submit_mpi_jobs
        submit_mpi_jobs_80x90
        submit_mpi_jobs_160x180
        echo "All jobs submitted successfully"
    else
        echo "Available options:"
        echo "1: Submit MPI jobs (40x40)"
        echo "2: Submit MPI jobs (80x90)"
        echo "3: Submit MPI jobs (160x180)"
        echo "Enter option number (or 'q' to quit):"
        
        read option
        case $option in
            #1) submit_mpi_jobs ;;
            2) submit_mpi_jobs_80x90 ;;
            3) submit_mpi_jobs_160x180 ;;
            q) echo "Exiting..." ;;
            *) echo "Invalid option" ;;
        esac
    fi
}

# 执行主函数
main