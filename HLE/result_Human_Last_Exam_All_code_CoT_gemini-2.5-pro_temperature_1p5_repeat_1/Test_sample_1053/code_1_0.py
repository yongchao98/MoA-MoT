import matplotlib.pyplot as plt

def solve_and_plot():
    """
    Calculates the path from a series of movements, plots the resulting shape,
    and prints the name of the shape.
    """
    # The list of relative movements (Δx, Δy)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initialize the starting point and the path
    x0, y0 = 0, 0
    path_points = [(x0, y0)]
    current_x, current_y = x0, y0

    # Calculate the path by iterating through the deltas
    for i, (dx, dy) in enumerate(deltas):
        next_x = current_x + dx
        next_y = current_y - dy # Note the subtraction of Δy

        # As requested, output the final equation numbers for the last point
        if i == len(deltas) - 1:
            print("Calculation for the final point:")
            print(f"x_final = {current_x:.2f} + {dx:.2f} = {next_x:.2f}")
            print(f"y_final = {current_y:.2f} - {dy:.2f} = {next_y:.2f}")
            print("-" * 20)

        path_points.append((next_x, next_y))
        current_x, current_y = next_x, next_y
        
    # The shape is identified from the plot
    shape_name = "cat"
    print(f"The shape is a: {shape_name}")

    # Close the shape for plotting by connecting the last point to the start
    plot_points = path_points + [(x0, y0)]

    # Unzip the points into separate x and y coordinate lists for plotting
    x_coords, y_coords = zip(*plot_points)

    # Plot the shape
    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')
    plt.title("Generated Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    # Use 'equal' aspect ratio to see the true shape
    plt.axis('equal')
    plt.grid(True)
    plt.show()

# Run the function
solve_and_plot()