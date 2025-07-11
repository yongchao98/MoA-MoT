import matplotlib.pyplot as plt

def solve_and_plot_path():
    """
    This function calculates and plots a path based on a series of relative movements,
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
    x, y = x0, y0
    
    # Store points for plotting
    points = [(x, y)]

    print("Calculating the path points:")
    # Calculate each subsequent point
    for dx, dy in deltas:
        prev_x, prev_y = x, y
        x = prev_x + dx
        y = prev_y - dy # Subtract dy to flip the y-axis
        points.append((x, y))
        # Output the calculation for the current step
        print(f"({prev_x:7.2f}, {prev_y:7.2f}) + ({dx:7.2f}, {-dy:7.2f}) = ({x:7.2f}, {y:7.2f})")

    # Print closing segment information
    last_x, last_y = points[-1]
    print("\nClosing the shape:")
    print(f"({last_x:7.2f}, {last_y:7.2f}) -> ({x0:7.2f}, {y0:7.2f})")
        
    # Close the shape by adding the starting point at the end for plotting
    points.append((x0, y0))
    
    # --- Plotting the shape ---
    # Unzip points into x and y coordinates
    x_coords, y_coords = zip(*points)

    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, 'b-o') # Blue line with circle markers
    
    # Mark the start point in green
    plt.plot(x0, y0, 'go', markersize=10, label='Start (0,0)')

    # Set plot properties
    plt.title('Path Plot')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box') # Ensure correct aspect ratio
    plt.legend()
    
    # To display the plot, uncomment the following line
    # plt.show()
    
    # Save the plot to a file to be able to see it without running plt.show()
    plt.savefig("shape_plot.png")
    print("\nPlot has been saved to shape_plot.png")

# Run the function
solve_and_plot_path()

# After observing the generated plot, the shape is a dog.
print("\nThe shape is a (single lowercase word): dog")