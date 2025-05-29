import matplotlib.pyplot as plt

# Define the points
points = [(41, 42), (37, 34), (41, 30), (45, 34), (41, 42)]

# Separate the x and y coordinates
x, y = zip(*points)

# Plot the shape
plt.figure()
plt.plot(x, y, marker='o')
plt.fill(x, y, alpha=0.3)
plt.title('SVG Path Shape')
plt.xlim(30, 50)
plt.ylim(25, 45)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()