import matplotlib.pyplot as plt

def plot_shape():
    """
    This function calculates and plots a shape based on a series of relative movements.
    """
    # The data for relative movements (dx, dy)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x0, y0 = 0.0, 0.0

    # Lists to store the coordinates of the path
    points = [(x0, y0)]
    current_x, current_y = x0, y0

    # Print the starting point and equation header
    print(f"Starting point (x_0, y_0) = ({x0:.2f}, {y0:.2f})")
    print("-" * 50)
    print("x_{i+1} = x_i + dx_i")
    print("y_{i+1} = y_i - dy_i")
    print("-" * 50)

    # Calculate the path by iterating through the deltas
    for i, (dx, dy) in enumerate(deltas):
        # Store previous point for printing
        prev_x, prev_y = current_x, current_y
        
        # Calculate the next point
        current_x = prev_x + dx
        current_y = prev_y - dy
        
        # Print the equation for the current step
        print(f"Step {i+1:2d}: (dx={dx:+.2f}, dy={dy:+.2f})")
        print(f"  x_{i+1} = {prev_x:.2f} + {dx:.2f} = {current_x:.2f}")
        print(f"  y_{i+1} = {prev_y:.2f} - {dy:.2f} = {current_y:.2f}")
        
        # Add the new point to our list
        points.append((current_x, current_y))

    # To close the shape, add the starting point at the end of the path
    points.append((x0, y0))
    print("\nClosing the shape by returning to the starting point (0.00, 0.00).")

    # Unzip the points into two separate lists for x and y coordinates
    x_coords, y_coords = zip(*points)

    # Plotting the shape
    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='purple')

    # Setting plot properties
    plt.title("The Plotted Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

# Run the function
plot_shape()

# By running the code, the plotted shape is clearly identifiable.
# The shape is a (cat).
