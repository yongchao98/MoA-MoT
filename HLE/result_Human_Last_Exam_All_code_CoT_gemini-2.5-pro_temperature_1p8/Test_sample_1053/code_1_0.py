import matplotlib.pyplot as plt

def solve_and_plot_path():
    """
    Calculates, prints, and plots a path from a series of movements,
    and identifies the resulting shape.
    """
    # The given relative movements (Δx, Δy)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x0, y0 = 0.0, 0.0
    
    # List to store all the points of the path
    points = [(x0, y0)]
    
    current_x, current_y = x0, y0
    
    print(f"Starting point (x_0, y_0) = ({current_x:.2f}, {current_y:.2f})\n")

    # Calculate the path point by point
    for i, (dx, dy) in enumerate(deltas):
        next_x = current_x + dx
        next_y = current_y - dy  # Note the subtraction of Δy

        print(f"Step {i+1}: Applying (Δx, Δy) = ({dx:.2f}, {dy:.2f})")
        print(f"  x_{i+1} = x_{i} + Δx_{i} = {current_x:.2f} + {dx:.2f} = {next_x:.2f}")
        print(f"  y_{i+1} = y_{i} - Δy_{i} = {current_y:.2f} - {dy:.2f} = {next_y:.2f}")
        print("-" * 30)

        points.append((next_x, next_y))
        current_x, current_y = next_x, next_y
    
    # Close the shape by connecting the last point to the starting point
    points.append((x0, y0))

    # Prepare data for plotting
    x_coords, y_coords = zip(*points)
    
    # Plot the final shape
    plt.figure(figsize=(7, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')
    plt.title("The Generated Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.axis('equal') # Ensure correct aspect ratio to identify the shape
    
    # Mark the start/end point in red
    plt.plot(x0, y0, 'ro', markersize=10, label='Start/End Point')
    plt.legend()
    
    plt.show()

# Run the function
solve_and_plot_path()
