import matplotlib.pyplot as plt

# The list of relative movements (Δx, Δy)
movements = [
    (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
    (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
    (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
    (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
    (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
]

# Initialize the starting point and the list of path points
x0, y0 = 0.0, 0.0
points = [(x0, y0)]
current_x, current_y = x0, y0

print("Calculating the path from the starting point (0.00, 0.00).")
print("The update rules are: x_new = x_old + Δx and y_new = y_old - Δy\n")

# Calculate the full path by iterating through the movements
for i, (dx, dy) in enumerate(movements):
    prev_x, prev_y = current_x, current_y
    current_x += dx
    current_y -= dy  # Subtract Δy to flip the y-axis
    points.append((current_x, current_y))
    
    # Output the calculation for the current step
    print(f"Step {i+1}: Applying movement (Δx={dx}, Δy={dy})")
    print(f"  x_{i+1} = {prev_x:.2f} + {dx} = {current_x:.2f}")
    print(f"  y_{i+1} = {prev_y:.2f} - {dy} = {current_y:.2f}")
    print(f"  New Point: ({current_x:.2f}, {current_y:.2f})\n")

# Close the shape by adding the starting point to the end of the list
points.append((x0, y0))
print("Path complete. Closing the shape by returning to the starting point (0.00, 0.00).")

# Prepare the x and y coordinates for plotting
x_coords = [p[0] for p in points]
y_coords = [p[1] for p in points]

# Plot the final shape
print("\nGenerating the plot...")
plt.figure(figsize=(8, 8))
plt.plot(x_coords, y_coords, marker='o', linestyle='-')
plt.title("Path Plot")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()