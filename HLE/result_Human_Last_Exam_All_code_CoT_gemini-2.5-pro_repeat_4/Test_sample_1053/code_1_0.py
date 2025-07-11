import matplotlib.pyplot as plt

def solve_and_plot():
    """
    Calculates the path points from a series of relative movements,
    prints the calculation for each step, and plots the resulting shape.
    """
    # The list of relative movements (delta_x, delta_y)
    delta_movements = [
        [0.15, 0.26], [0.38, 0.10], [1.32, -0.06], [0.53, -0.07],
        [1.06, -1.02], [0, 1.32], [-1.06, -1.06], [-0.26, -0.13],
        [-1.12, -0.05], [0.28, -0.48], [0.02, -0.15], [-0.21, 0.26],
        [-0.03, -0.34], [-0.14, 0.26], [-0.17, -0.19], [0.23, 0.38],
        [0.02, 0.25], [-0.78, 0.17], [-0.22, 0.26], [0, 0.26]
    ]

    # Initial starting point
    x0, y0 = 0.0, 0.0
    
    # List to store the points of the path
    points = [(x0, y0)]
    x_current, y_current = x0, y0

    print(f"Starting at (x_0, y_0) = ({x_current:.2f}, {y_current:.2f})\n")

    # Calculate the path points
    for i, (dx, dy) in enumerate(delta_movements):
        x_next = x_current + dx
        # Note the subtraction of dy to flip the y-axis
        y_next = y_current - dy
        
        print(f"Step {i+1}:")
        print(f"  Current point (x_{i}, y_{i}) = ({x_current:.2f}, {y_current:.2f})")
        print(f"  Movement (Δx_{i}, Δy_{i}) = ({dx}, {dy})")
        print(f"  x_{i+1} = {x_current:.2f} + {dx} = {x_next:.2f}")
        print(f"  y_{i+1} = {y_current:.2f} - {dy} = {y_next:.2f}")
        print(f"  New point (x_{i+1}, y_{i+1}) = ({x_next:.2f}, {y_next:.2f})\n")

        points.append((x_next, y_next))
        x_current, y_current = x_next, y_next

    # Close the shape by adding the starting point to the end of the list
    points.append((x0, y0))

    # Unzip the points into x and y coordinate lists for plotting
    x_coords, y_coords = zip(*points)

    # Plot the shape
    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-')
    plt.title("Generated Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    # Ensure the aspect ratio is equal to avoid distortion
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Save the plot to a file
    plot_filename = "shape.png"
    plt.savefig(plot_filename)
    
    print(f"The plot has been saved as '{plot_filename}'.")
    print("The shape is a cat.")

if __name__ == '__main__':
    solve_and_plot()