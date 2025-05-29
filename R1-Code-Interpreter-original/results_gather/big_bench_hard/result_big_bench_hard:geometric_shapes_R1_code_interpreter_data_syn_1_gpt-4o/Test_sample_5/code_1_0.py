import matplotlib.pyplot as plt

# Define the points
points = [
    (2.10, 61.72),
    (38.86, 44.98),
    (65.79, 63.32),
    (39.63, 11.62),
    (21.25, 37.04),
    (31.13, 38.33),
    (2.10, 61.72)  # Closing the path
]

# Plot the points and lines
plt.figure()
for i in range(len(points) - 1):
    x_values = [points[i][0], points[i+1][0]]
    y_values = [points[i][1], points[i+1][1]]
    plt.plot(x_values, y_values, 'bo-')

# Set the aspect of the plot to be equal
plt.gca().set_aspect('equal', adjustable='box')
plt.title('SVG Path Visualization')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()