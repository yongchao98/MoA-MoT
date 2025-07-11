import matplotlib.pyplot as plt

def plot_shape():
    """
    Calculates and plots a path based on a series of relative movements,
    then identifies the resulting shape.
    """
    # The data for relative movements (delta_x, delta_y)
    movements = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x0, y0 = 0, 0
    x_coords = [x0]
    y_coords = [y0]

    # Current point, initialized to the starting point
    x_i, y_i = x0, y0

    print(f"Starting point (x_0, y_0) = ({x0:.2f}, {y0:.2f})\n")

    # Calculate the path and print each step's equation
    for i, (delta_x, delta_y) in enumerate(movements):
        print(f"Step {i+1}: Calculating point (x_{i+1}, y_{i+1})")
        
        x_iplus1 = x_i + delta_x
        # Note the subtraction of delta_y to flip the y-axis
        y_iplus1 = y_i - delta_y

        # Print the equations for the current step
        print(f"x_{i+1} = x_{i} + \u0394x_{i} = {x_i:.2f} + {delta_x:.2f} = {x_iplus1:.2f}")
        print(f"y_{i+1} = y_{i} - \u0394y_{i} = {y_i:.2f} - {delta_y:.2f} = {y_iplus1:.2f}")
        print(f"New point (x_{i+1}, y_{i+1}) = ({x_iplus1:.2f}, {y_iplus1:.2f})\n")

        x_coords.append(x_iplus1)
        y_coords.append(y_iplus1)
        
        # Update the current point for the next iteration
        x_i, y_i = x_iplus1, y_iplus1

    # Close the shape by adding the starting point at the end
    x_coords.append(x0)
    y_coords.append(y0)
    print("Closing shape by drawing a line segment back to the starting point (0.00, 0.00).")

    # Create the plot
    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-')

    # Set plot properties
    plt.title("Path Plot")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    # Ensure correct aspect ratio to avoid distortion
    plt.gca().set_aspect('equal', adjustable='box')

    # Show the plot
    plt.show()

# Run the function
plot_shape()