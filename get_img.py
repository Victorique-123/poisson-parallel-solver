import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_poisson(input_file, M, N):
    plt.rcParams['font.family'] = 'DejaVu Sans'

    # Read data from TSV file
    data = np.loadtxt(input_file, delimiter='\t')
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # Get grid dimensions
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    X, Y = np.meshgrid(x_unique, y_unique)
    Z = z.reshape(len(y_unique), len(x_unique))

    # Create figure with subplots
    fig = plt.figure(figsize=(15, 6))

    # 3D Surface plot
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(X, Y, Z, 
                           cmap='coolwarm',
                           linewidth=0.2,
                           antialiased=True)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('u(x,y)')
    ax1.set_title(f'Трехмерный график решения ({M}×{N})')
    fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)

    # Set 3D view angle
    ax1.view_init(elev=30, azim=45)

    # 2D Contour plot
    ax2 = fig.add_subplot(122)
    contour = ax2.contourf(X, Y, Z, 
                          levels=20, 
                          cmap='coolwarm')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title(f'Контурный график решения ({M}×{N})')
    fig.colorbar(contour, ax=ax2)

    # Add grid points
    ax2.plot(X, Y, 'k.', markersize=1, alpha=0.3)

    # Add title
    plt.suptitle(f'Численное решение уравнения Пуассона на сетке {M}×{N}', y=1.05)

    # Adjust layout and save
    plt.tight_layout()
    output_file = f'3D_solution_{M}x{N}.png'
    plt.savefig(output_file, 
                dpi=300, 
                bbox_inches='tight')
    print(f"The image has been saved as: {output_file}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='绘制泊松方程数值解的可视化图形')
    parser.add_argument('input_file', help='输入的TSV文件路径')
    parser.add_argument('M', type=int, help='网格M的大小')
    parser.add_argument('N', type=int, help='网格N的大小')
    
    args = parser.parse_args()
    
    plot_poisson(args.input_file, args.M, args.N)

if __name__ == "__main__":
    main()
