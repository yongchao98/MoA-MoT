import matplotlib.pyplot as plt
import numpy as np

def solve_and_plot():
    """
    Calculates the path coordinates from relative movements,
    prints the coordinates, and plots the resulting shape.
    """
    # The table of relative movements (delta_x, delta_y)
    movements = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x0, y0 = 0.0, 0.0
    points = [(x0, y0)]

    # Calculate the path by iterating through movements
    current_x, current_y = x0, y0
    for dx, dy in movements:
        next_x = current_x + dx
        # Subtract dy to flip the y-axis
        next_y = current_y - dy
        points.append((next_x, next_y))
        current_x, current_y = next_x, next_y

    # Close the shape by adding the starting point to the end of the path
    points.append((x0, y0))

    # Print the coordinates of each point in the path
    print("The coordinates of the path vertices are:")
    for i, (x, y) in enumerate(points):
        # Using np.round to handle potential floating point inaccuracies for cleaner output
        print(f"({np.round(x, 2):.2f}, {np.round(y, 2):.2f})")

    # Prepare data for plotting
    x_coords, y_coords = zip(*points)

    # Plot the shape
    plt.figure(figsize=(7, 7))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')
    plt.title("Shape from Path Segments")
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    # Ensure the aspect ratio is equal to avoid distortion
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()

solve_and_plot()