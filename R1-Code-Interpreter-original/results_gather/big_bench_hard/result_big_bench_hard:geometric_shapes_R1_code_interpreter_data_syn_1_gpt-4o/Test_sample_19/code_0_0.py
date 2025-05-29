import matplotlib.pyplot as plt

# Define the points
points = [
    (64.00, 63.00),
    (44.00, 63.00),
    (44.00, 50.00),
    (64.00, 50.00),
    (64.00, 45.00),
    (85.00, 57.00),
    (64.00, 68.00),
    (64.00, 63.00)  # Closing the path
]

# Separate the x and y coordinates
x, y = zip(*points)

# Plot the shape
plt.figure()
plt.plot(x, y, marker='o')
plt.fill(x, y, alpha=0.3)
plt.title('SVG Path Shape')
plt.xlim(40, 90)
plt.ylim(40, 70)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()