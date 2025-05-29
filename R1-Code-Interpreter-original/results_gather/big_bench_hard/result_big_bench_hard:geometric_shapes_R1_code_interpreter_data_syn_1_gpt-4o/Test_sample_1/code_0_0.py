import matplotlib.pyplot as plt

# Define the points
points = [(30.17, 45.97), (58.79, 40.36), (18.10, 15.70), (30.17, 45.97)]

# Unzip the points into x and y coordinates
x, y = zip(*points)

# Plot the points and lines
plt.figure()
plt.plot(x, y, marker='o')
plt.title('SVG Path Visualization')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()