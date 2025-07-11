import matplotlib.pyplot as plt

def solve_and_plot_path():
    """
    Calculates a path from a series of relative movements, plots it,
    and prints the full calculation for the final point.
    """
    # The initial starting point
    x0, y0 = 0.0, 0.0

    # The series of relative movements (delta x, delta y)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initialize path list and current coordinates
    points = [(x0, y0)]
    current_x, current_y = x0, y0
    
    # Initialize strings to build the full equations for the final point
    x_eq_str = f"{x0:.2f}"
    y_eq_str = f"{y0:.2f}"

    # Calculate the path by iterating through the deltas
    for dx, dy in deltas:
        # Update current coordinates
        current_x += dx
        current_y -= dy  # Note the subtraction of delta_y to flip the y-axis
        points.append((current_x, current_y))
        
        # Append the current movement to the x-equation string
        if dx >= 0:
            x_eq_str += f" + {dx:.2f}"
        else:
            x_eq_str += f" - {abs(dx):.2f}"
            
        # Append the current movement to the y-equation string
        # The y-update is y = y - dy, so we are adding the (-dy) term
        term_y = -dy
        if term_y >= 0:
            y_eq_str += f" + {term_y:.2f}"
        else:
            y_eq_str += f" - {abs(term_y):.2f}"

    # The last calculated point is the final point of the open path
    final_x = current_x
    final_y = current_y

    # Print the full equations showing how each number is used
    print("Final X coordinate calculation:")
    print(f"x_final = {x_eq_str} = {final_x:.2f}")
    print("\nFinal Y coordinate calculation:")
    print(f"y_final = {y_eq_str} = {final_y:.2f}")
    
    # Close the shape by adding the starting point to the end of the path
    points.append((x0, y0))

    # Unzip the list of point tuples into separate lists for x and y coordinates
    x_coords, y_coords = zip(*points)

    # Plot the final closed shape
    plt.figure(figsize=(7, 7))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')
    plt.title("Plot of the Path")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.grid(True)
    # Use 'equal' aspect ratio for accurate shape representation
    plt.axis('equal')
    plt.show()

# Run the function
solve_and_plot_path()