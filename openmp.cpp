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
#include <omp.h>
#include <algorithm> 

using namespace std;

// 求解器性能统计结构体
struct SolverStats {
    size_t iterations;
    double solve_time;
    double delta_achieved;
};

// 使用一维向量实现的网格数据结构
struct Grid {
    vector<double> values;
    size_t global_size_x, global_size_y;
    double step_x, step_y, a1, b1, a2, b2;
    
    // 访问器方法
    double& at(size_t i, size_t j) {
        return values[i * (global_size_y + 1) + j];
    }
    
    const double& at(size_t i, size_t j) const {
        return values[i * (global_size_y + 1) + j];
    }
};

// 创建并初始化网格
Grid create_grid(double a1_, double a2_, double b1_, double b2_, size_t m, size_t n) {
    Grid grid;
    grid.global_size_x = m;
    grid.global_size_y = n;
    grid.a1 = a1_;
    grid.b1 = b1_;
    grid.a2 = a2_;
    grid.b2 = b2_;
    grid.step_x = (b1_ - a1_) / static_cast<double>(m);
    grid.step_y = (b2_ - a2_) / static_cast<double>(n);
    grid.values.resize((m + 1) * (n + 1), 0.0);
    return grid;
}

// 坐标计算函数
double get_x(const Grid& grid, size_t i) {
    return grid.a1 + i * grid.step_x;
}

double get_y(const Grid& grid, size_t j) {
    return grid.a2 + j * grid.step_y;
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

// 计算系数aij和bij
void calculate_aij(Grid& result, double eps) {
    double inv_h2 = 1.0 / result.step_y;
    
    #pragma omp parallel for collapse(2) schedule(guided)
    for (size_t i = 1; i <= result.global_size_x; ++i) {
        for (size_t j = 1; j <= result.global_size_y; ++j) {
            double x = get_x(result, i);
            double y = get_y(result, j);
            double l = get_y_length(x - 0.5 * result.step_x, 
                                  y - 0.5 * result.step_y, 
                                  y + 0.5 * result.step_y);
            result.at(i, j) = inv_h2 * l + (1.0 - inv_h2 * l) / eps;
        }
    }
}

void calculate_bij(Grid& result, double eps) {
    double inv_h1 = 1.0 / result.step_x;
    
    #pragma omp parallel for collapse(2) schedule(guided)
    for (size_t i = 1; i <= result.global_size_x; ++i) {
        for (size_t j = 1; j <= result.global_size_y; ++j) {
            double x = get_x(result, i);
            double y = get_y(result, j);
            double l = get_x_length(y - 0.5 * result.step_y,
                                  x - 0.5 * result.step_x, 
                                  x + 0.5 * result.step_x);
            result.at(i, j) = inv_h1 * l + (1.0 - inv_h1 * l) / eps;
        }
    }
}

// 应用差分算子
void apply_diff_operator(Grid& result, const Grid& gf, const Grid& aij, const Grid& bij) {
    #pragma omp parallel for collapse(2) schedule(guided)
    for (size_t i = 1; i < gf.global_size_x; ++i) {
        for (size_t j = 1; j < gf.global_size_y; ++j) {
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
}

// 计算标量积
double scalar_product(const Grid& gf1, const Grid& gf2) {
    double result = 0.0;
    double weight = gf1.step_x * gf1.step_y;
    
    #pragma omp parallel for collapse(2) reduction(+:result) schedule(guided)
    for (size_t i = 0; i < gf1.global_size_x; ++i) {
        for (size_t j = 0; j < gf1.global_size_y; ++j) {
            result += weight * gf1.at(i, j) * gf2.at(i, j);
        }
    }
    return result;
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
    // 找到基准点（最左下角的点）
    Point ref = findReferencePoint(points);
    
    // 将基准点移到第一个位置
    for (int i = 0; i < points.size(); i++) {
        if (points[i].x == ref.x && points[i].y == ref.y) {
            swap(points[i], points[0]);
            break;
        }
    }
    
    // 按极角排序其余的点
    sort(points.begin() + 1, points.end(),
        [ref](const Point& a, const Point& b) {
            double angleA = polarAngle(ref, a);
            double angleB = polarAngle(ref, b);
            if (angleA != angleB) return angleA < angleB;
            // 如果极角相同，距离近的点排在前面
            return distance(ref, a) < distance(ref, b);
        });
}

// 计算多边形面积
double calculatePolygonArea(vector<Point>& points) {
    int n = points.size();
    if (n < 3) return 0;  // 少于3个点无法构成多边形
    
    // 先对点进行极角排序
    sortPoints(points);
    
    double area = 0.0;
    // 使用向量叉积计算面积
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += points[i].x * points[j].y;
        area -= points[j].x * points[i].y;
    }
    
    return fabs(area) / 2.0;
}



double calculate_area(double x1, double x2, double y1, double y2) {
    //首先处理与直线x=1的关系
    if (x1 >= 1.0) return 0.0;
    //限制右边界不超过x=1
    double x2_new = min(x2, 1.0);

    // 计算抛物线与网格单元的交点
    vector<Point> points;
    
    // 添加网格单元的四个顶点(如果在区域内)
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
    
    // 添加抛物线与网格单元边界的交点
    // 与垂直边界的交点
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
        
    // 与水平边界的交点
    double x_intersect = y1*y1;
    if (x_intersect >= x1 && x_intersect <= x2_new) {
        points.push_back(Point(x_intersect, y1));
    }
    x_intersect = y2*y2;
    if (x_intersect >= x1 && x_intersect <= x2_new) {
        points.push_back(Point(x_intersect, y2));
    }
    
    // 如果没有交点,说明完全在区域外
    if (points.empty()) return 0.0;
    
    // 使用多边形面积计算函数计算面积
    return calculatePolygonArea(points);
}



// 修改后的init_grid函数
void init_grid(Grid& gf) {
    #pragma omp parallel for collapse(2) schedule(guided)
    for (size_t i = 1; i < gf.global_size_x; ++i) {
        for (size_t j = 1; j < gf.global_size_y; ++j) {
            double x1 = get_x(gf, i) - 0.5 * gf.step_x;
            double x2 = get_x(gf, i) + 0.5 * gf.step_x;
            double y1 = get_y(gf, j) - 0.5 * gf.step_y;
            double y2 = get_y(gf, j) + 0.5 * gf.step_y;
            
            double area = calculate_area(x1, x2, y1, y2);
            gf.at(i, j) = area / ((x2-x1)*(y2-y1));
        }
    }
    // 边界点设为0
    for (size_t i = 0; i <= gf.global_size_x; ++i) {
        gf.at(i, 0) = 0.0;
        gf.at(i, gf.global_size_y) = 0.0;
    }
    for (size_t j = 0; j <= gf.global_size_y; ++j) {
        gf.at(0, j) = 0.0;
        gf.at(gf.global_size_x, j) = 0.0;
    }
}

// 主求解函数
SolverStats solve_poisson(double a1, double a2, double b1, double b2,
                         size_t m, size_t n, double delta_stop,
                         const string& output_file, int num_threads) {
    SolverStats stats = {0, 0.0, 0.0};
    double total_start = omp_get_wtime();
    
    omp_set_num_threads(num_threads);
    
    Grid w = create_grid(a1, a2, b1, b2, m, n);
    Grid r = create_grid(a1, a2, b1, b2, m, n);  // 也是 Aw
    Grid ar = create_grid(a1, a2, b1, b2, m, n);
    Grid b = create_grid(a1, a2, b1, b2, m, n);
    Grid aij = create_grid(a1, a2, b1, b2, m, n);
    Grid bij = create_grid(a1, a2, b1, b2, m, n);


    double h_max = max(w.step_x, w.step_y);
    double epsilon = max(h_max * h_max, 1.00e-06);
    
    init_grid(b);
    calculate_aij(aij, epsilon);
    calculate_bij(bij, epsilon);

    double solve_start = omp_get_wtime();
    
    int iteration_count = 0;
    while (true) {
        iteration_count++;
        apply_diff_operator(r, w, aij, bij);
        
        #pragma omp parallel for collapse(2) schedule(guided)
        for (size_t i = 0; i < w.global_size_x; ++i) {
            for (size_t j = 0; j < w.global_size_y; ++j) {
                r.at(i, j) = r.at(i, j) - b.at(i, j);
            }
        }
        
        double r_norm = scalar_product(r, r);
        apply_diff_operator(ar, r, aij, bij);
        double ar_r = scalar_product(ar, r);
        double tau = r_norm / (ar_r + 1e-10);


        //计算解的增量的能量范数
        double delta_w = 0.0;
        #pragma omp parallel for collapse(2) reduction(+:delta_w) schedule(guided)
        for (size_t i = 0; i < w.global_size_x; ++i) {
            for (size_t j = 0; j < w.global_size_y; ++j) {
                double diff = tau * r.at(i, j);
                delta_w += diff * diff * w.step_x * w.step_y;
            }
        }
        stats.delta_achieved = sqrt(delta_w); 
        // stats.delta_achieved = sqrt(r_norm) * fabs(tau);       
        if (stats.delta_achieved <= delta_stop) break;
        
        
        #pragma omp parallel for collapse(2) schedule(guided)
        for (size_t i = 0; i < w.global_size_x; ++i) {
            for (size_t j = 0; j < w.global_size_y; ++j) {
                w.at(i, j) -= tau * r.at(i, j);
            }
        }

        
        ++stats.iterations;
    }
    //储存计算时间
    stats.solve_time = omp_get_wtime() - solve_start;
    
    // if (!output_file.empty()) {
    //     #pragma omp master
    //     {
    //         ofstream fout(output_file.c_str());
    //         if (!fout) {
    //             throw runtime_error("Failed to open output file");
    //         } else {
    //             for (size_t i = 0; i <= w.size_m; ++i) {
    //                 for (size_t j = 0; j <= w.size_n; ++j) {
    //                     double x = get_x(w, i);
    //                     double y = get_y(w, j);
    //                     fout << x << "\t" << y << "\t" << w.at(i,j) << "\n";
    //                 }
    //             }
    //         }
    //     }
    // }
    return stats;
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

int string_to_int(const char* str) {
    stringstream ss(str);
    int result;
    ss >> result;
    if (ss.fail()) {
        throw runtime_error("Invalid number format");
    }
    return result;
}

int main(int argc, char **argv) {
    try {
        if (argc != 5) {
            cerr << "Usage: " << argv[0] << " m n delta_w_stop num_threads" << endl;
            return 1;
        }

        size_t M = string_to_size_t(argv[1]);
        size_t N = string_to_size_t(argv[2]);
        double delta_w_stop = string_to_double(argv[3]);
        int num_threads = string_to_int(argv[4]);

        stringstream ss;
        ss << "opm_" << M << "x" << N <<".tsv";
        string output_file = ss.str();

        cout << "Solving Poisson equation (OpenMP)" << endl;
        cout << "Grid Size: " << M << "x" << N << endl;
        cout << "Number of Threads: " << num_threads << endl;
        cout << "Convergence Delta: " << delta_w_stop << endl;

        SolverStats stats = solve_poisson(0.0, -1.0, 1.0, 1.0, M, N, 
                                        delta_w_stop, output_file, num_threads);

        cout << "\nResults:" << endl;
        cout << "Iterations: " << stats.iterations << endl;
        cout << "Final delta: " << scientific << setprecision(3) << stats.delta_achieved << endl;
        cout << "Solve time: " << fixed << setprecision(3) << stats.solve_time << " seconds" << endl;

        ofstream time_file("opm_times.tsv", ios::app);
        if (time_file) {
            time_file << "M:" << M << " | "
                    << "N:" << N << " | " 
                    << "num_threads:" << num_threads << " | "
                    << "delta_w_stop:" << delta_w_stop << " | "
                    << "Final delta:" << scientific << setprecision(3) << stats.delta_achieved << " | "
                    << "Iterations:" << stats.iterations << " | "
                    << "Solve time:" << fixed << setprecision(3) << stats.solve_time << "\n";
        }

        return 0;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}