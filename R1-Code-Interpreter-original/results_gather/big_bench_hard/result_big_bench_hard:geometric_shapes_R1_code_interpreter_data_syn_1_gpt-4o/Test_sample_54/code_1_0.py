import matplotlib.pyplot as plt

# Define the points
points = [
    (64.00, 63.00),
    (44.00, 63.00),
    (44.00, 50.00),
    (64.00, 50.00),
    (64.00, 37.00),
    (85.00, 57.00),
    (64.00, 76.00),
    (64.00, 63.00)  # Closing the shape
]

# Separate the x and y coordinates
x, y = zip(*points)

# Plot the shape
plt.figure()
plt.plot(x, y, marker='o')
plt.fill(x, y, alpha=0.3)
plt.title('SVG Path Shape')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.axis('equal')
plt.show()