import matplotlib.pyplot as plt

def solve_and_plot():
    """
    Calculates the path from a series of movements, prints each step,
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
    
    print("Calculating path points:")
    # Calculate each subsequent point
    for i, (dx, dy) in enumerate(deltas):
        x_prev, y_prev = x, y
        x = x_prev + dx
        y = y_prev - dy
        points.append((x, y))
        print(f"Step {i+1}: ({x:.2f}, {y:.2f}) = ({x_prev:.2f}, {y_prev:.2f}) + ({dx:.2f}, -{dy:.2f})")

    # Final step to close the shape by returning to the origin
    x_prev, y_prev = x, y
    x_start, y_start = points[0]
    dx_close = x_start - x_prev
    dy_close = y_prev - y_start
    points.append((x_start, y_start))
    print(f"Closing Step: ({x_start:.2f}, {y_start:.2f}) = ({x_prev:.2f}, {y_prev:.2f}) + ({dx_close:.2f}, -{dy_close:.2f})")


    # Prepare data for plotting
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]

    # Plot the shape
    plt.figure(figsize=(6, 8))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-')
    plt.title("The Generated Shape")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    # Ensure the aspect ratio is equal to see the true shape
    plt.gca().set_aspect('equal', adjustable='box')
    print("\nDisplaying plot. Close the plot window to continue.")
    plt.show()

    print("\nThe shape is a cat.")

if __name__ == '__main__':
    solve_and_plot()