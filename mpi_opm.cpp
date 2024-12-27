#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cassert>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <mpi.h>
#include <omp.h>

using namespace std;

struct SolverStats {
    size_t iterations;
    double solve_time;
    double delta_achieved;
};
//ceshi
// 扩展Grid结构以支持MPI
struct Grid {
    vector<double> values;
    size_t size_m, size_n;           // 全局大小
    size_t local_m, local_n;         // 局部大小
    size_t offset_m, offset_n;       // 在全局网格中的偏移
    double step_x, step_y, a1, b1, a2, b2;
    
    int coords[2];                   // 进程在网格中的坐标
    int dims[2];                     // 进程网格的维度
    int neighbors[4];                // 相邻进程的rank (left, right, top, bottom)
    MPI_Comm cart_comm;             // 笛卡尔通信器
    
    vector<double> send_buf[4];      // 发送缓冲区
    vector<double> recv_buf[4];      // 接收缓冲区
    
    double& at(size_t i, size_t j) {
        return values[i * (local_n + 2) + j];
    }
    
    const double& at(size_t i, size_t j) const {
        return values[i * (local_n + 2) + j];
    }
};

// 首先添加新的辅助函数在 create_grid_mpi 函数之前
void calculate_local_sizes(size_t total_size, int num_procs, int proc_idx, 
                         size_t& local_size, size_t& offset) {
    size_t base_size = total_size / num_procs;
    size_t remainder = total_size % num_procs;
    
    // 如果有余数，前remainder个进程多分配一个节点
    if (proc_idx < remainder) {
        local_size = base_size + 1;
        offset = proc_idx * (base_size + 1);
    } else {
        local_size = base_size;
        offset = remainder * (base_size + 1) + (proc_idx - remainder) * base_size;
    }
}

Grid create_grid_mpi(double a1_, double a2_, double b1_, double b2_, 
                    size_t m, size_t n, int size, int rank) {
    Grid grid;
    
    // 设置全局参数
    grid.size_m = m;
    grid.size_n = n;
    grid.a1 = a1_;
    grid.b1 = b1_;
    grid.a2 = a2_;
    grid.b2 = b2_;
    grid.step_x = (b1_ - a1_) / m;
    grid.step_y = (b2_ - a2_) / n;
    
    // 创建笛卡尔拓扑
    int dims[2] = {0, 0};
    int periods[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    
    // 保存进程网格维度
    grid.dims[0] = dims[0];
    grid.dims[1] = dims[1];
    
    // 创建笛卡尔通信器
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid.cart_comm);
    MPI_Cart_coords(grid.cart_comm, rank, 2, grid.coords);
    
    // 使用新的划分方法计算局部大小和偏移
    calculate_local_sizes(m, dims[0], grid.coords[0], grid.local_m, grid.offset_m);
    calculate_local_sizes(n, dims[1], grid.coords[1], grid.local_n, grid.offset_n);
    
    // 分配内存（包括ghost cells）
    grid.values.resize((grid.local_m + 2) * (grid.local_n + 2), 0.0);
    
    // 获取相邻进程
    MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.neighbors[0], &grid.neighbors[1]);
    MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.neighbors[2], &grid.neighbors[3]);
    
    // 初始化通信缓冲区
    for (int i = 0; i < 4; ++i) {
        grid.send_buf[i].resize(i < 2 ? grid.local_n : grid.local_m);
        grid.recv_buf[i].resize(i < 2 ? grid.local_n : grid.local_m);
    }
    
    return grid;
}


// 交换边界数据
void exchange_boundaries(Grid& grid) {
    MPI_Request requests[8];
    MPI_Status statuses[8];
    int req_count = 0;
    
    // 发送/接收左右边界
    if(grid.neighbors[0] != MPI_PROC_NULL) {
        for(size_t j = 0; j < grid.local_n; ++j) {
            grid.send_buf[0][j] = grid.at(1, j+1);
        }
        MPI_Isend(grid.send_buf[0].data(), grid.local_n, MPI_DOUBLE, 
                  grid.neighbors[0], 0, grid.cart_comm, &requests[req_count++]);
        MPI_Irecv(grid.recv_buf[0].data(), grid.local_n, MPI_DOUBLE,
                  grid.neighbors[0], 1, grid.cart_comm, &requests[req_count++]);
    }
    
    if(grid.neighbors[1] != MPI_PROC_NULL) {
        for(size_t j = 0; j < grid.local_n; ++j) {
            grid.send_buf[1][j] = grid.at(grid.local_m, j+1);
        }
        MPI_Isend(grid.send_buf[1].data(), grid.local_n, MPI_DOUBLE,
                  grid.neighbors[1], 1, grid.cart_comm, &requests[req_count++]);
        MPI_Irecv(grid.recv_buf[1].data(), grid.local_n, MPI_DOUBLE,
                  grid.neighbors[1], 0, grid.cart_comm, &requests[req_count++]);
    }
    
    // 发送/接收上下边界
    if(grid.neighbors[2] != MPI_PROC_NULL) {
        for(size_t i = 0; i < grid.local_m; ++i) {
            grid.send_buf[2][i] = grid.at(i+1, 1);
        }
        MPI_Isend(grid.send_buf[2].data(), grid.local_m, MPI_DOUBLE,
                  grid.neighbors[2], 2, grid.cart_comm, &requests[req_count++]);
        MPI_Irecv(grid.recv_buf[2].data(), grid.local_m, MPI_DOUBLE,
                  grid.neighbors[2], 3, grid.cart_comm, &requests[req_count++]);
    }
    
    if(grid.neighbors[3] != MPI_PROC_NULL) {
        for(size_t i = 0; i < grid.local_m; ++i) {
            grid.send_buf[3][i] = grid.at(i+1, grid.local_n);
        }
        MPI_Isend(grid.send_buf[3].data(), grid.local_m, MPI_DOUBLE,
                  grid.neighbors[3], 3, grid.cart_comm, &requests[req_count++]);
        MPI_Irecv(grid.recv_buf[3].data(), grid.local_m, MPI_DOUBLE,
                  grid.neighbors[3], 2, grid.cart_comm, &requests[req_count++]);
    }
    
    // 等待所有通信完成
    MPI_Waitall(req_count, requests, statuses);
    
    // 更新ghost cells
    if(grid.neighbors[0] != MPI_PROC_NULL) {
        for(size_t j = 0; j < grid.local_n; ++j) {
            grid.at(0, j+1) = grid.recv_buf[0][j];
        }
    }
    
    if(grid.neighbors[1] != MPI_PROC_NULL) {
        for(size_t j = 0; j < grid.local_n; ++j) {
            grid.at(grid.local_m+1, j+1) = grid.recv_buf[1][j];
        }
    }
    
    if(grid.neighbors[2] != MPI_PROC_NULL) {
        for(size_t i = 0; i < grid.local_m; ++i) {
            grid.at(i+1, 0) = grid.recv_buf[2][i];
        }
    }
    
    if(grid.neighbors[3] != MPI_PROC_NULL) {
        for(size_t i = 0; i < grid.local_m; ++i) {
            grid.at(i+1, grid.local_n+1) = grid.recv_buf[3][i];
        }
    }
}

// 并行版本的标量积计算
double scalar_product(const Grid& gf1, const Grid& gf2) {
    double local_result = 0.0;
    double global_result = 0.0;
    double weight = gf1.step_x * gf1.step_y;
    #pragma omp parallel for reduction(+:local_result) collapse(2)
    for (size_t i = 1; i <= gf1.local_m; ++i) {
        for (size_t j = 1; j <= gf1.local_n; ++j) {
            local_result += weight * gf1.at(i, j) * gf2.at(i, j);
        }
    }

    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, gf1.cart_comm);

    return global_result;
}

// 辅助函数实现
double get_x_local(const Grid& grid, size_t i) {
    return grid.a1 + (grid.offset_m + i) * grid.step_x;
}

double get_y_local(const Grid& grid, size_t j) {
    return grid.a2 + (grid.offset_n + j) * grid.step_y;
}

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

// 计算两点之间的距离
double distance(const Point& p1, const Point& p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

// 计算向量的叉积
double crossProduct(const Point& p1, const Point& p2, const Point& p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

// 有限差分导数计算
double der_x_right(const Grid& grid, size_t i, size_t j) {
    return (grid.at(i+1, j) - grid.at(i, j)) / grid.step_x;
}

double der_x_left(const Grid& grid, size_t i, size_t j) {
    return (grid.at(i, j) - grid.at(i-1, j)) / grid.step_x;
}

double der_y_right(const Grid& grid, size_t i, size_t j) {
    return (grid.at(i, j+1) - grid.at(i, j)) / grid.step_y;
}

double der_y_left(const Grid& grid, size_t i, size_t j) {
    return (grid.at(i, j) - grid.at(i, j-1)) / grid.step_y;
}

double get_parabolic_y(double x) {
    return (x > 0.0) ? sqrt(x) : -1.0;
}

// 计算y方向长度
double get_y_length(double x, double y1, double y2) {
    if (x >= 1.0 || x <= 0.0) return 0.0;
    double y_parabola = get_parabolic_y(x);
    if (y_parabola < 0.0) return 0.0;
    
    double y_min = -y_parabola;
    double y_max = y_parabola;
    
    if (y1 > y_max || y2 < y_min) return 0.0;
    return min(y2, y_max) - max(y1, y_min);
}

double get_x_length(double y, double x1, double x2) {
    double x_min = y * y;
    double x_max = 1.0;
    
    if (x1 > x_max || x2 < x_min) return 0.0;
    return min(x2, x_max) - max(x1, x_min);
}

// 找到最左下角的点（作为基准点）
Point findReferencePoint(vector<Point>& points) {
    Point ref = points[0];
    for (const Point& p : points) {
        if (p.y < ref.y || (p.y == ref.y && p.x < ref.x)) {
            ref = p;
        }
    }
    return ref;
}

// 计算点相对于基准点的极角
double polarAngle(const Point& ref, const Point& p) {
    double dx = p.x - ref.x;
    double dy = p.y - ref.y;
    return atan2(dy, dx);
}

// 按极角排序点
void sortPoints(vector<Point>& points) {
    Point ref = findReferencePoint(points);
    
    for (int i = 0; i < points.size(); i++) {
        if (points[i].x == ref.x && points[i].y == ref.y) {
            swap(points[i], points[0]);
            break;
        }
    }
    
    sort(points.begin() + 1, points.end(),
        [ref](const Point& a, const Point& b) {
            double angleA = polarAngle(ref, a);
            double angleB = polarAngle(ref, b);
            if (angleA != angleB) return angleA < angleB;
            return distance(ref, a) < distance(ref, b);
        });
}

// 计算多边形面积
double calculatePolygonArea(vector<Point>& points) {
    int n = points.size();
    if (n < 3) return 0;
    
    sortPoints(points);
    
    double area = 0.0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += points[i].x * points[j].y;
        area -= points[j].x * points[i].y;
    }
    
    return fabs(area) / 2.0;
}


double calculate_area(double x1, double x2, double y1, double y2) {
    if (x1 >= 1.0) return 0.0;
    double x2_new = min(x2, 1.0);

    vector<Point> points;
    
    if (y1*y1 < x1) {
        points.push_back(Point(x1, y1));
    }
    if (y1*y1 < x2_new) {
        points.push_back(Point(x2_new, y1));
    }
    if (y2*y2 < x1) {
        points.push_back(Point(x1, y2));
    }
    if (y2*y2 < x2_new) {
        points.push_back(Point(x2_new, y2));
    }
    
    double y_intersect;
    y_intersect = sqrt(x1);
    if (y_intersect >= y1 && y_intersect <= y2) {
        points.push_back(Point(x1, y_intersect));
    }
    y_intersect = -sqrt(x1);
    if (y_intersect >= y1 && y_intersect <= y2) {
        points.push_back(Point(x1, -y_intersect));
    }
    y_intersect = sqrt(x2_new);
    if (y_intersect >= y1 && y_intersect <= y2) {
        points.push_back(Point(x2_new, y_intersect));
    }
    y_intersect = -sqrt(x2_new);
    if (y_intersect >= y1 && y_intersect <= y2) {
        points.push_back(Point(x2_new, -y_intersect));
    }
        
    double x_intersect = y1*y1;
    if (x_intersect >= x1 && x_intersect <= x2_new) {
        points.push_back(Point(x_intersect, y1));
    }
    x_intersect = y2*y2;
    if (x_intersect >= x1 && x_intersect <= x2_new) {
        points.push_back(Point(x_intersect, y2));
    }
    
    if (points.empty()) return 0.0;
    
    return calculatePolygonArea(points);
}


void init_grid_local(Grid& gf) {
    // 首先计算内部节点的值
    #pragma omp parallel for collapse(2)
    for (size_t i =0; i < gf.local_m; ++i) {
        for (size_t j = 0; j < gf.local_n; ++j) {
            double x = get_x_local(gf, i);
            double y = get_y_local(gf, j);
            
            double x1 = x - 0.5 * gf.step_x;
            double x2 = x + 0.5 * gf.step_x;
            double y1 = y - 0.5 * gf.step_y;
            double y2 = y + 0.5 * gf.step_y;
            double area = calculate_area(x1, x2, y1, y2);
            gf.at(i+1, j+1) = area / ((x2-x1)*(y2-y1));
        }
    }

    // 检查是否是最左边的进程(coords[0] == 0)
    if (gf.coords[0] == 0) {
        // 设置左边第二列（索引为1）为0
        for (size_t j = 1; j <= gf.local_n; ++j) {
            gf.at(1, j) = 0.0;
        }
    }

    // 检查是否是最上边的进程(coords[1] == 0)
    if (gf.coords[1] == 0) {
        // 设置上边第二行（索引为1）为0
        for (size_t i = 1; i <= gf.local_m; ++i) {
            gf.at(i, 1) = 0.0;
        }
    }
}

void calculate_aij_local(Grid& result, double eps) {
    double inv_h2 = 1.0 / result.step_y;
    
    #pragma omp parallel for collapse(2)
    for(size_t i = 0; i <= result.local_m; ++i) {
        for(size_t j = 0; j <= result.local_n; ++j) {
            double x = get_x_local(result, i);
            double y = get_y_local(result, j);
            double l = get_y_length(x - 0.5 * result.step_x,
                                  y - 0.5 * result.step_y,
                                  y + 0.5 * result.step_y);
            result.at(i+1, j+1) = inv_h2 * l + (1.0 - inv_h2 * l) / eps;
        }
    }
}

void calculate_bij_local(Grid& result, double eps) {
    double inv_h1 = 1.0 / result.step_x;
    
    #pragma omp parallel for collapse(2)
    for(size_t i = 0; i <= result.local_m; ++i) {
        for(size_t j = 0; j <= result.local_n; ++j) {
            double x = get_x_local(result, i);
            double y = get_y_local(result, j);
            double l = get_x_length(y - 0.5 * result.step_y,
                                  x - 0.5 * result.step_x,
                                  x + 0.5 * result.step_x);
            result.at(i+1, j+1) = inv_h1 * l + (1.0 - inv_h1 * l) / eps;
        }
    }
}

void apply_diff_operator_local(Grid& result, const Grid& gf,
                             const Grid& aij, const Grid& bij) {

    #pragma omp parallel for collapse(2)
    for(size_t i = 1; i <= gf.local_m; ++i) {
        for(size_t j = 1; j <= gf.local_n; ++j) {
            double a2 = aij.at(i + 1, j);
            double a1 = aij.at(i, j);
            double b2 = bij.at(i, j + 1);
            double b1 = bij.at(i, j);

            result.at(i, j) = -(
                (a2 * der_x_right(gf, i, j) - a1 * der_x_left(gf, i, j)) / gf.step_x +
                (b2 * der_y_right(gf, i, j) - b1 * der_y_left(gf, i, j)) / gf.step_y
            );
        }
    }
    // 检查是否是最左边的进程(coords[0] == 0)
    if (result.coords[0] == 0) {
        // 设置左边第二列（索引为1）为0
        for (size_t j = 0; j <= result.local_n + 1; ++j) {
            result.at(1, j) = 0.0;
        }
    }

    // 检查是否是最上边的进程(coords[1] == 0)
    if (result.coords[1] == 0) {
        // 设置上边第二行（索引为1）为0
        for (size_t i = 0; i <= result.local_m + 1; ++i) {
            result.at(i, 1) = 0.0;
        }
    }
}

// 字符串转换工具
size_t string_to_size_t(const char* str) {
    stringstream ss(str);
    size_t result;
    ss >> result;
    if (ss.fail()) {
        throw runtime_error("Invalid number format");
    }
    return result;
}

double string_to_double(const char* str) {
    stringstream ss(str);
    double result;
    ss >> result;
    if (ss.fail()) {
        throw runtime_error("Invalid number format");
    }
    return result;
}


SolverStats solve_poisson_mpi(double a1, double a2, double b1, double b2,
                             size_t m, size_t n, double delta_stop) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    SolverStats stats = {0, 0.0, 0.0};
    double start_time = MPI_Wtime();
    
    Grid w = create_grid_mpi(a1, a2, b1, b2, m, n, size, rank);
    Grid r = create_grid_mpi(a1, a2, b1, b2, m, n, size, rank);
    Grid ar = create_grid_mpi(a1, a2, b1, b2, m, n, size, rank);
    Grid b = create_grid_mpi(a1, a2, b1, b2, m, n, size, rank);
    Grid aij = create_grid_mpi(a1, a2, b1, b2, m, n, size, rank);
    Grid bij = create_grid_mpi(a1, a2, b1, b2, m, n, size, rank);
    
    double h_max = max(w.step_x, w.step_y);
    double epsilon = max(h_max * h_max, 1.0e-6);
    
    init_grid_local(b);
    calculate_aij_local(aij, epsilon);
    calculate_bij_local(bij, epsilon);

    int iteration_count = 0;
    while (true) {
        iteration_count++;
        
        exchange_boundaries(w);

        apply_diff_operator_local(r, w, aij, bij);
        
        #pragma omp parallel for collapse(2)
        for (size_t i = 1; i <= w.local_m; ++i) {
            for (size_t j = 1; j <= w.local_n; ++j) {
                r.at(i, j) = r.at(i, j) - b.at(i, j);
            }
        }
        
        double r_norm = scalar_product(r, r);
        
        exchange_boundaries(r);
        apply_diff_operator_local(ar, r, aij, bij);
        
        double ar_r = scalar_product(ar, r);
        double tau = r_norm / (ar_r + 1e-10);
        
        double local_delta_w = 0.0;
        #pragma omp parallel for reduction(+:local_delta_w) collapse(2)
        for (size_t i = 1; i <= w.local_m; ++i) {
            for (size_t j = 1; j <= w.local_n; ++j) {
                double diff = tau * r.at(i, j);
                local_delta_w += diff * diff * w.step_x * w.step_y;
            }
        }
        
        double global_delta_w;
        MPI_Allreduce(&local_delta_w, &global_delta_w, 1, MPI_DOUBLE, MPI_SUM, w.cart_comm);
        stats.delta_achieved = sqrt(global_delta_w);
        
        if (stats.delta_achieved <= delta_stop) {
            break;
        }
        
        #pragma omp parallel for collapse(2)
        for (size_t i = 1; i <= w.local_m; ++i) {
            for (size_t j = 1; j <= w.local_n; ++j) {
                w.at(i, j) -= tau * r.at(i, j);
            }
        }
        
        ++stats.iterations;
    }
    
    stats.solve_time = MPI_Wtime() - start_time;
    return stats;
}



int main(int argc, char **argv) {
    try {
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
        
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        
        if(argc != 5 && rank == 0) {
            cerr << "Usage: " << argv[0] << " m n delta_w_stop thread" << endl;
            MPI_Finalize();
            return 1;
        }

        size_t M = string_to_size_t(argv[1]);
        size_t N = string_to_size_t(argv[2]);
        double delta_w_stop = string_to_double(argv[3]);
        int num_thread = string_to_size_t(argv[4]);

        omp_set_num_threads(num_thread);

        // stringstream ss;
        // ss << "mpi_" << M << "x" << N << "_p" << size << ".tsv";
        // string output_file = ss.str();

        if(rank == 0) {
            cout << "Solving Poisson equation (MPI)" << endl;
            cout << "Grid Size: " << M << "x" << N << endl;
            cout << "Number of Processes: " << size << endl;
            cout << "Number of Threads: " <<  num_thread << endl;
            cout << "Convergence Delta: " << delta_w_stop << endl;
        }

        SolverStats stats = solve_poisson_mpi(0.0, -1.0, 1.0, 1.0, M, N,
                                            delta_w_stop);

        if(rank == 0) {
            cout << "\nResults:" << endl;
            cout << "Iterations: " << stats.iterations << endl;
            cout << "Final delta: " << scientific << setprecision(3) 
                 << stats.delta_achieved << endl;
            cout << "Solve time: " << fixed << setprecision(3) 
                 << stats.solve_time << " seconds" << endl;

            ofstream time_file("mpi_times.tsv", ios::app);
            if(time_file) {
                time_file << "M:" << M << " | "
                         << "N:" << N << " | "
                         << "Processes:" << size << " | "
                         << "Threads:" << num_thread << " | "
                         << "delta_w_stop:" << delta_w_stop << " | "
                         << "Final delta:" << scientific << setprecision(3) 
                         << stats.delta_achieved << " | "
                         << "Iterations:" << stats.iterations << " | "
                         << "Solve time:" << fixed << setprecision(3) 
                         << stats.solve_time << "\n";
            }
        }

        MPI_Finalize();
        return 0;
    }
    catch(const exception& e) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) {
            cerr << "Error: " << e.what() << endl;
        }
        MPI_Finalize();;
        return 1;
    }
}