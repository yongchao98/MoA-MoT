import matplotlib.pyplot as plt

def solve_and_plot():
    """
    This function calculates and plots a path based on a series of movements,
    and identifies the resulting shape.
    """
    # Define the relative movements (Δx, Δy)
    deltas = [
        [0.15, 0.26], [0.38, 0.10], [1.32, -0.06], [0.53, -0.07],
        [1.06, -1.02], [0, 1.32], [-1.06, -1.06], [-0.26, -0.13],
        [-1.12, -0.05], [0.28, -0.48], [0.02, -0.15], [-0.21, 0.26],
        [-0.03, -0.34], [-0.14, 0.26], [-0.17, -0.19], [0.23, 0.38],
        [0.02, 0.25], [-0.78, 0.17], [-0.22, 0.26], [0, 0.26]
    ]

    # Initial starting point
    x0, y0 = 0, 0
    points = [(x0, y0)]
    x_current, y_current = x0, y0

    print(f"Starting at point P0 = ({x_current:.2f}, {y_current:.2f})\n")
    print("Calculating path segments:")

    # Calculate the path by iterating through the deltas
    for i, (dx, dy) in enumerate(deltas):
        x_prev, y_prev = x_current, y_current
        
        # Compute the next point
        x_current = x_prev + dx
        y_current = y_prev - dy  # Note the subtraction of dy
        
        points.append((x_current, y_current))
        
        # Print the calculation for the current segment
        print(f"P{i} -> P{i+1}: x = {x_prev:.2f} + {dx:.2f} = {x_current:.2f}, y = {y_prev:.2f} - {dy:.2f} = {y_current:.2f}")

    # Add the starting point at the end to close the shape
    print("\nClosing the shape by connecting the last point back to the start.")
    print(f"P{len(deltas)} -> P0: From ({x_current:.2f}, {y_current:.2f}) to ({x0:.2f}, {y0:.2f})")
    points.append((x0, y0))

    # Prepare data for plotting
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]

    # Plot the shape
    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-')
    plt.title("Path Plot")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.axis('equal') # Use an equal aspect ratio
    plt.show()

# Run the function
solve_and_plot()