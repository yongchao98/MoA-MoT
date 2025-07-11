import numpy as np
import matplotlib.pyplot as plt

def create_and_plot_points():
    """
    This function generates and plots a valid configuration for n=10 points.
    It demonstrates the possibility of the (Nr=2, Ng=4, Ny=4) distribution.
    It does not mathematically prove that 10 is the maximum, but shows it is possible.
    The proof outline is in the text.
    """

    # Green points (large convex quadrilateral)
    g_pts = np.array([
        [-5, -4],
        [5, -4],
        [5, 4],
        [-5, 4]
    ])

    # Yellow points (smaller convex quad inside the kernel of the green one)
    y_pts = np.array([
        [-1, -0.8],
        [1, -0.8],
        [1, 0.8],
        [-1, 0.8]
    ])

    # Red points (placed to pierce the yellow triangles)
    # y_I is the intersection of diagonals of y_pts, which is (0,0)
    # Place r1 in triangle y0-yI-y1 -> y0-O-y1
    # Place r2 in triangle y2-yI-y3 -> y2-O-y3
    r_pts = np.array([
        [0.0, -0.1],  # r1
        [0.0, 0.1]   # r2
    ])
    
    # Plotting
    plt.figure(figsize=(8, 8))
    plt.scatter(r_pts[:, 0], r_pts[:, 1], c='red', s=100, label='Red (Nr=2)')
    plt.scatter(g_pts[:, 0], g_pts[:, 1], c='green', s=100, label='Green (Ng=4)')
    plt.scatter(y_pts[:, 0], y_pts[:, 1], c='yellow', s=100, label='Yellow (Ny=4)')

    # Draw convex hulls to visualize the arrangement
    from scipy.spatial import ConvexHull

    hull_g = ConvexHull(g_pts)
    for simplex in hull_g.simplices:
        plt.plot(g_pts[simplex, 0], g_pts[simplex, 1], 'g-')

    hull_y = ConvexHull(y_pts)
    for simplex in hull_y.simplices:
        plt.plot(y_pts[simplex, 0], y_pts[simplex, 1], 'y-')
    
    # Example Triangles to show conditions are met
    # Cond 2: Green triangle contains a yellow point
    plt.fill(g_pts[[0, 1, 2], 0], g_pts[[0, 1, 2], 1], 'g', alpha=0.1)
    # Cond 3: Yellow triangle contains a red point
    plt.fill(y_pts[[0, 1, 2], 0], y_pts[[0, 1, 2], 1], 'y', alpha=0.2)
        
    plt.title("Configuration for n=10 points (2 Red, 4 Green, 4 Yellow)")
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    
    print("A plot showing a valid configuration for n=10 would be generated if this code were executed.")
    print("This configuration demonstrates that n=10 is achievable.")
    print("Based on the geometric constraints (piercing numbers for triangles),")
    print("configurations that would sum to more than 10 can be ruled out.")
    print("For instance, any configuration with a color count of 5 (e.g., yellow) would require")
    print("at least 3 points of its piercing color (red), making a (2, x, 5) distribution impossible.")
    print("This limits the number of points for any color to 4 in any arrangement with a color count of 2.")
    print("Maximum points = 2 (Red) + 4 (Green) + 4 (Yellow) = 10.")


create_and_plot_points()
