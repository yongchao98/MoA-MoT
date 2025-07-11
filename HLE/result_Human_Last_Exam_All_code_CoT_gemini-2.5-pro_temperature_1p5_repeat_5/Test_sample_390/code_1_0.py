import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Let's consider the case n=3, d=3 for visualization.
# The vectors {y_i} are linearly independent. We will demonstrate the shape
# for the fundamental case where {y_i} form an orthonormal basis, as this
# reveals the underlying geometric structure.
# Let y1, y2, y3 be the standard basis vectors.
y1 = np.array([1., 0., 0.])
y2 = np.array([0., 1., 0.])
y3 = np.array([0., 0., 1.])
basis_vectors = [y1, y2, y3]

# The vector 's' is any unit vector in the span of {y_i}, which is R^3 in this case.
# We will generate a large number of random unit vectors on the surface of a sphere.
num_points = 5000
# Generate points from a 3D Gaussian distribution and normalize them to get
# a uniform distribution on the sphere.
s_vectors = np.random.normal(size=(num_points, 3))
s_vectors /= np.linalg.norm(s_vectors, axis=1)[:, np.newaxis]

# For each vector 's', we compute the corresponding point in the set S.
# A point in S is (|<y1,s>|^2, |<y2,s>|^2, |<y3,s>|^2).
S_points = np.zeros((num_points, 3))
for i, y in enumerate(basis_vectors):
    S_points[:, i] = np.abs(np.dot(s_vectors, y))**2

# For an orthonormal basis, the sum of the components of any point in S
# should equal ||s||^2, which is 1. This is a key property of a simplex.
# Let's verify this and print the result for a few points.
sums_of_components = np.sum(S_points, axis=1)

print("Let's test the case where the basis {y_i} is orthonormal.")
print("The points in S are (|<y1,s>|^2, |<y2,s>|^2, |<y3,s>|^2).")
print("Their components must sum to 1.\n")
for i in range(5):
    p = S_points[i]
    s = sums_of_components[i]
    # The output shows each component of the point and their sum.
    print(f"Point {i+1}: ({p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f}). Sum = {s:.6f}")

print("\nAs the sums are consistently 1 (within floating point precision), all points")
print("lie on the plane x1 + x2 + x3 = 1.")
print("Since each component x_i >= 0, the set S is precisely the standard 2-simplex.")

# Visualization of the simplex
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
# Plot the points of S
ax.scatter(S_points[:, 0], S_points[:, 1], S_points[:, 2], s=5, alpha=0.5, label='Points in S')

# Define the vertices of the 2-simplex (triangle)
vertices = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
# Create a Poly3DCollection for the simplex
simplex = Poly3DCollection([vertices], alpha=0.2, facecolor='cyan', edgecolor='k')
ax.add_collection3d(simplex)

# Set plot labels and title
ax.set_title("The Shape of Set S for n=3 (Orthonormal Case)")
ax.set_xlabel("$x_1 = |\langle y_1, s \\rangle|^2$")
ax.set_ylabel("$x_2 = |\langle y_2, s \\rangle|^2$")
ax.set_zlabel("$x_3 = |\langle y_3, s \\rangle|^2$")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_zlim(0, 1)
# Set aspect ratio to be equal
ax.set_box_aspect([1, 1, 1])
ax.view_init(elev=30, azim=45)
plt.show()
