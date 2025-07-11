import matplotlib.pyplot as plt

def generate_and_plot_shape():
    """
    Calculates path coordinates from relative movements, prints them,
    and plots the resulting shape.
    """
    # Table of relative movements (Δx, Δy)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x, y = 0.0, 0.0
    points = [(x, y)]

    print(f"Starting Point 0: ({x:.2f}, {y:.2f})")

    # Calculate the path
    for i, (dx, dy) in enumerate(deltas):
        x += dx
        # Subtract dy to flip the y-axis
        y -= dy
        points.append((x, y))
        print(f"Point {i+1}: ({x:.2f}, {y:.2f})")

    # Close the shape by connecting back to the start
    points.append(points[0])

    # Unzip points for plotting
    x_coords, y_coords = zip(*points)

    # Plot the shape
    plt.figure(figsize=(6, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-')
    plt.title("Path Plot")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.axis('equal') # Ensure aspect ratio is 1:1
    plt.show()

# Execute the function
generate_and_plot_shape()

print("\nThe shape is a cat.")