import matplotlib.pyplot as plt

# Define the points
points = [
    (4.70, 85.15),
    (2.36, 37.14),
    (18.76, 36.34),
    (21.10, 84.35),
    (4.70, 85.15)  # Closing the path
]

# Unzip the points into x and y coordinates
x, y = zip(*points)

# Plot the points
plt.figure()
plt.plot(x, y, marker='o')
plt.title('SVG Path Visualization')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()