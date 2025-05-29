import matplotlib.pyplot as plt

# Define the points
points = [(41, 42), (37, 32), (41, 30), (45, 32), (41, 42)]

# Separate the x and y coordinates
x, y = zip(*points)

# Plot the shape
plt.figure()
plt.plot(x, y, marker='o')
plt.title('SVG Path Shape')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.axis('equal')  # Equal scaling for x and y axes
plt.show()