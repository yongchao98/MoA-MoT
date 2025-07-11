import matplotlib.pyplot as plt

def solve_and_plot_path():
    """
    Calculates and plots a path from a series of movements,
    and prints the calculation for each step.
    """
    # The list of relative movements (dx, dy)
    movements = [
        [0.15, 0.26], [0.38, 0.10], [1.32, -0.06], [0.53, -0.07],
        [1.06, -1.02], [0, 1.32], [-1.06, -1.06], [-0.26, -0.13],
        [-1.12, -0.05], [0.28, -0.48], [0.02, -0.15], [-0.21, 0.26],
        [-0.03, -0.34], [-0.14, 0.26], [-0.17, -0.19], [0.23, 0.38],
        [0.02, 0.25], [-0.78, 0.17], [-0.22, 0.26], [0, 0.26]
    ]

    # Initial starting point
    x0, y0 = 0.0, 0.0

    # Lists to store the coordinates of the path
    x_coords = [x0]
    y_coords = [y0]

    # Current position, initialized to the starting point
    xi, yi = x0, y0

    print("--- Path Calculation Steps ---")
    # Calculate the path by iterating through the movements
    for i, (dx, dy) in enumerate(movements):
        x_prev, y_prev = xi, yi
        # Compute the next point
        xi += dx
        yi -= dy
        
        # Print the equation for the current step
        print(f"x_{i+1} = {x_prev:.2f} + {dx:.2f} = {xi:.2f}")
        print(f"y_{i+1} = {y_prev:.2f} - {dy:.2f} = {yi:.2f}")

        # Store the new point
        x_coords.append(xi)
        y_coords.append(yi)

    # Close the shape by adding a line back to the starting point
    x_coords.append(x0)
    y_coords.append(y0)

    print("\n--- Plotting Shape ---")
    print("The shape drawn by the path is a cat.")

    # Plot the path
    plt.figure(figsize=(7, 7))
    plt.plot(x_coords, y_coords, 'b-') # 'b-' creates a blue line
    plt.plot(x_coords, y_coords, 'ro', markersize=4) # 'ro' creates red dots for points
    plt.title("Path that forms a shape")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.grid(True)
    # Ensures the aspect ratio is 1:1, crucial for shape recognition
    plt.axis('equal')
    plt.show()

# Run the function to perform calculations and generate the plot
solve_and_plot_path()