import matplotlib.pyplot as plt

def solve_and_plot_shape():
    """
    Calculates the path coordinates, plots the resulting shape,
    and prints the name of the shape.
    """
    # The relative movements (dx, dy)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initialize the starting point and the list of path points
    x0, y0 = 0, 0
    points = [(x0, y0)]
    current_x, current_y = x0, y0

    # Calculate the coordinates of the path
    for dx, dy in deltas:
        next_x = current_x + dx
        # Subtract dy to flip the y-axis as per the instruction
        next_y = current_y - dy
        points.append((next_x, next_y))
        current_x, current_y = next_x, next_y

    # Close the shape by returning to the starting point
    points.append((x0, y0))

    # Identify the shape
    shape_name = "bird"
    print(f"The shape is a {shape_name}.")

    # Prepare data for plotting
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]

    # Plot the shape
    plt.figure(figsize=(8, 8))
    # Plot line segments
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')
    # Highlight the starting point
    plt.plot(x0, y0, marker='*', markersize=15, color='r', label='Start/End Point')
    
    plt.title("Plot of the Generated Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    # Use 'equal' aspect ratio to see the true shape
    plt.axis('equal')
    plt.legend()
    plt.show()

# Run the function
solve_and_plot_shape()