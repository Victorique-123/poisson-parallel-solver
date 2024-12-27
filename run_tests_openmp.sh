#!/bin/bash


# 加载必要的模块
module load SpectrumMPI/10.1.0

# 编译所有程序
compile_programs() {
    echo "Compiling programs..."
    g++ -fopenmp openmp.cpp -o openmp_O3 -std=c++11 -O3
    
    if [ $? -eq 0 ]; then
        echo "Compilation successful"
    else
        echo "Compilation failed"
        exit 1
    fi
}


# 提交OpenMP程序作业（O3优化版本）
submit_openmp_O3_jobs() {
    echo "Submitting OpenMP O3 jobs for 40x40 grid..."
    
    for threads in 1 2 4 8 16; do
        bsub -n 1 -q normal \
             -o opm_O3_40x40.log \
             -e error.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J openmp_O3_${threads}t \
             -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
             "OMP_NUM_THREADS=4 ./openmp_O3 40 40 1.00e-07 ${threads}"
    done
}

# 提交OpenMP CG作业 (80x90)
submit_openmp_O3_jobs_80x90() {
    echo "Submitting OpenMP CG jobs for 80x90 grid..."
    
    for threads in 1 2 4 8 16 32; do
        bsub -n 1 -q normal \
             -o opm_80x90.log \
             -e error_opm_80x90.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J openmp_O3_80x90_${threads}t \
             -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
             "OMP_NUM_THREADS=4 ./openmp_O3 80 90 1.00e-08 ${threads}"
    done
}

# 提交OpenMP CG作业 (160x180)
submit_openmp_O3_jobs_160x180() {
    echo "Submitting OpenMP CG jobs for 160x180 grid..."
    
    for threads in 4 8 16 32 1 2; do
        bsub -n 1 -q normal \
             -o opm_160x180.log \
             -e error_opm_160x180.log \
             -m "polus-c3-ib polus-c4-ib" \
             -J openmp_cg_160x180_${threads}t \
             -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
             "OMP_NUM_THREADS=4 ./openmp_O3 160 180 2.50e-09 ${threads}"
    done
}

# 主函数
main() {
    compile_programs
    
    echo "Do you want to submit all jobs? [y/n]"
    read answer
    
    if [ "$answer" = "y" ]; then
        submit_openmp_O3_jobs
        submit_openmp_O3_jobs_80x90
        submit_openmp_O3_jobs_160x180
        echo "All jobs submitted successfully"
    else
        echo "Available options:"
        echo "4: Submit OpenMP O3 jobs (40x40)"
        echo "5: Submit OpenMP O3 jobs (80x90)"
        echo "6: Submit OpenMP O3 jobs (160x180)"
        echo "Enter option number (or 'q' to quit):"
        
        read option
        case $option in
            1) submit_openmp_O3_jobs ;;
            2) submit_openmp_O3_jobs_80x90 ;;
            3) submit_openmp_O3_jobs_160x180 ;;
            q) echo "Exiting..." ;;
            *) echo "Invalid option" ;;
        esac
    fi
}

# 执行主函数
main