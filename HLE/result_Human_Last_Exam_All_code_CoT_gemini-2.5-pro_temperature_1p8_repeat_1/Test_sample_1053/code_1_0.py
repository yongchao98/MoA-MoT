import matplotlib.pyplot as plt

# The list of relative movements (Δx, Δy)
deltas = [
    (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
    (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
    (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
    (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
    (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
]

# Initial starting point (x₀, y₀)
x_points = [0]
y_points = [0]

# Current point, initialized to the starting point
x_i, y_i = 0, 0

print(f"Starting at (x_0, y_0) = ({x_i:.2f}, {y_i:.2f})\n")

# Calculate each subsequent point
for i, (dx, dy) in enumerate(deltas):
    x_next = x_i + dx
    y_next = y_i - dy
    
    # Output the calculation for each point as per the formula
    print(f"x_{i+1} = {x_i:.2f} + {dx:.2f} = {x_next:.2f}")
    print(f"y_{i+1} = {y_i:.2f} - {dy:.2f} = {y_next:.2f}")
    print("-" * 20)
    
    x_points.append(x_next)
    y_points.append(y_next)
    
    # Update the current point for the next iteration
    x_i, y_i = x_next, y_next

# Close the shape by returning to the starting point
x_points.append(x_points[0])
y_points.append(y_points[0])

# Plotting the shape
plt.figure(figsize=(8, 8))
plt.plot(x_points, y_points, marker='o', linestyle='-')
plt.title("Path Plot")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.grid(True)
plt.axis('equal') # Ensure aspect ratio is 1:1 to avoid distortion
plt.show()

# Based on the plot, the shape is a pig.
print("The shape is a pig.")