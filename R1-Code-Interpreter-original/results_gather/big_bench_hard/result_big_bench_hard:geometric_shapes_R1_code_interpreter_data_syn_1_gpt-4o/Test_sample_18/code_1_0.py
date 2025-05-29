import matplotlib.pyplot as plt

# Define the points
points = [
    (63.28, 25.83),
    (51.43, 14.65),
    (46.46, 19.91),
    (58.31, 31.09),
    (63.28, 25.83)
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