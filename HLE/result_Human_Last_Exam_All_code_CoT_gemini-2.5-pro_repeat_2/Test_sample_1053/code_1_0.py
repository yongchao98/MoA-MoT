import matplotlib.pyplot as plt

def solve_and_plot():
    """
    Calculates the path from the given movements, identifies the shape,
    and prints the required information.
    """
    # The list of relative movements (delta_x, delta_y)
    movements = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initialize starting point
    x0, y0 = 0, 0
    points = [(x0, y0)]
    
    x_current, y_current = x0, y0
    
    # Calculate all subsequent points
    for dx, dy in movements:
        x_current += dx
        y_current -= dy  # Subtract dy to flip the y-axis
        points.append((round(x_current, 2), round(y_current, 2)))
        
    # The last point of the path before closing
    last_point = points[-1]
    
    # The identified shape
    shape_name = "cat"
    print(f"The shape is a {shape_name}.")

    # The "final equation" is interpreted as the closing line segment
    # connecting the last point back to the start.
    print(f"The final closing segment connects the point ({last_point[0]}, {last_point[1]}) back to the starting point ({x0}, {y0}).")
    
    # The following code is for visualization and can be uncommented to see the plot
    # To close the shape for plotting, add the starting point to the end
    # points.append((x0, y0))
    # x_coords, y_coords = zip(*points)
    # plt.figure(figsize=(6,6))
    # plt.plot(x_coords, y_coords, 'b-')
    # plt.title("The Identified Shape")
    # plt.xlabel("X")
    # plt.ylabel("Y")
    # plt.grid(True)
    # plt.axis('equal')
    # plt.show()

solve_and_plot()