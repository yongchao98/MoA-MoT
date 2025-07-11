import numpy as np
import matplotlib.pyplot as plt

def get_shape_points(y_vectors, num_points=400):
    """
    Calculates the points in the set S for a given set of y_vectors.
    s is parameterized on the unit circle in the span of y_vectors.
    For n=d=2, the span is R^2 and s can be (cos(t), sin(t)).
    """
    thetas = np.linspace(0, 2 * np.pi, num_points)
    s_vectors = np.array([np.cos(thetas), np.sin(thetas)]).T
    
    points = []
    for s in s_vectors:
        z = tuple(np.abs(np.dot(y, s))**2 for y in y_vectors)
        points.append(z)
        
    return np.array(points)

# Case 1: Orthogonal vectors
y_ortho = [np.array([1.5, 0]), np.array([0, 1])]
points_ortho = get_shape_points(y_ortho)
# The equation for this case is z_1/1.5^2 + z_2/1^2 = 1 => z_1/2.25 + z_2 = 1
# This is a line segment (a 1-simplex).

# Case 2: Non-orthogonal vectors
y_non_ortho = [np.array([1, 0]), np.array([1, 1])]
points_non_ortho = get_shape_points(y_non_ortho)
# We derived the equation for this case as 4*z_1^2 + z_2^2 - 4*z_1 - 2*z_2 + 1 = 0.
# This is an ellipse.

# Print the final equations for clarity
print("Equation for the orthogonal case: z_1 / (1.5*1.5) + z_2 / (1*1) = 1.0")
print("Equation for the non-orthogonal case: 4*z_1^2 + z_2^2 - 4*z_1 - 2*z_2 + 1 = 0")


# Plotting the results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot for orthogonal case
ax1.plot(points_ortho[:, 0], points_ortho[:, 1], 'r-')
ax1.set_title('Shape S for Orthogonal Vectors (A Simplex)')
ax1.set_xlabel('$z_1 = |\langle y_1, s \rangle|^2$')
ax1.set_ylabel('$z_2 = |\langle y_2, s \rangle|^2$')
ax1.grid(True)
ax1.axis('equal')

# Plot for non-orthogonal case
ax2.plot(points_non_ortho[:, 0], points_non_ortho[:, 1], 'b-')
ax2.set_title('Shape S for Non-Orthogonal Vectors (An Ellipse)')
ax2.set_xlabel('$z_1 = |\langle y_1, s \rangle|^2$')
ax2.set_ylabel('$z_2 = |\langle y_2, s \rangle|^2$')
ax2.grid(True)
ax2.axis('equal')

plt.tight_layout()
plt.show()
