import matplotlib.pyplot as plt

# Define the points
points = [(41.97, 21.34), (12.23, 35.25), (4.62, 18.98), (34.36, 5.07), (41.97, 21.34)]

# Separate the x and y coordinates
x_coords, y_coords = zip(*points)

# Plot the points and lines
plt.figure()
plt.plot(x_coords, y_coords, marker='o')
plt.title('SVG Path Visualization')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()