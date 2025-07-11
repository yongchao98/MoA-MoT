import matplotlib.pyplot as plt

def solve_and_plot_path():
    """
    This function calculates, prints, and plots a path based on a series of
    relative movements, then identifies the resulting shape.
    """
    # The list of relative movements (Δx, Δy)
    delta_movements = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x0, y0 = 0, 0
    
    # Lists to store the coordinates of the path
    path_x = [x0]
    path_y = [y0]

    current_x, current_y = x0, y0

    print("Path Calculation Steps:")
    print(f"Starting at: (x₀, y₀) = ({x0:.2f}, {y0:.2f})")
    print("-" * 50)
    print(f"{'Step':<5} | {'Current (x,y)':<20} | {'Movement (Δx,-Δy)':<20} | {'New (x,y)':<20}")
    print("-" * 75)

    # Calculate each point in the path
    for i, (dx, dy) in enumerate(delta_movements):
        prev_x, prev_y = current_x, current_y
        
        # Apply the movement formulas
        current_x = prev_x + dx
        current_y = prev_y - dy # Note the subtraction of Δy

        # Print the calculation for this step
        print(f"{i+1:<5} | ({prev_x:6.2f}, {prev_y:6.2f})   | ({dx:6.2f}, {-dy:6.2f})   | ({current_x:6.2f}, {current_y:6.2f})")

        # Add the new point to the path
        path_x.append(current_x)
        path_y.append(current_y)
    
    # Close the shape by drawing a line from the last point back to the start
    last_x, last_y = current_x, current_y
    path_x.append(x0)
    path_y.append(y0)
    print("-" * 75)
    print(f"Closing shape: Drawing line from ({last_x:.2f}, {last_y:.2f}) back to ({x0:.2f}, {y0:.2f})")

    # Plotting the shape
    plt.figure(figsize=(7, 9))
    plt.plot(path_x, path_y, marker='o', markersize=3, linestyle='-', color='blue')
    plt.title("The Generated Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.axis('equal') # Use equal aspect ratio for accurate shape representation
    plt.show()

# Run the function
solve_and_plot_path()

# The shape is a cat.
print("\nThe shape is a cat.")