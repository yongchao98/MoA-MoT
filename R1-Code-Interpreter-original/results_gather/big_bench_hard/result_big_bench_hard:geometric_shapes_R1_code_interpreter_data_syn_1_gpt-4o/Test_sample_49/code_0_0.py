import matplotlib.pyplot as plt

# Define the points
points = [
    (81.97, 51.68),
    (56.45, 29.26),
    (50.05, 29.03),
    (49.41, 34.81),
    (54.75, 40.20),
    (81.97, 51.68)
]

# Plot the points and lines
plt.figure()
for i in range(len(points) - 1):
    x_values = [points[i][0], points[i+1][0]]
    y_values = [points[i][1], points[i+1][1]]
    plt.plot(x_values, y_values, marker='o')

# Set plot limits
plt.xlim(45, 85)
plt.ylim(25, 55)

# Show the plot
plt.title("SVG Path Visualization")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()