import numpy as np

def linf_norm(v):
    """Calculates the L-infinity norm of a vector."""
    return np.max(np.abs(v))

# We use the example of R^3 with the L-infinity norm.
# The unit ball for this norm is known to be metrically convex.
# The unit sphere is the surface of the cube [-1, 1]^3.

# We define two points p and q on an edge of the unit sphere.
p = np.array([1, 1, 1])
q = np.array([-1, 1, 1])

# Verify that p and q are on the unit sphere.
norm_p = linf_norm(p)
norm_q = linf_norm(q)
print(f"Point p = {p}, Norm ||p||_inf = {norm_p}")
print(f"Point q = {q}, Norm ||q||_inf = {norm_q}")
print("Both points lie on the unit sphere.\n")

# A point on the line segment is x(t) = (1-t)p + tq = (1-2t, 1, 1).
# The norm of any point on the segment is ||x(t)||_inf = max(|1-2t|, 1, 1) = 1.
# This confirms the entire segment lies on the unit sphere.

# Calculate the length of the segment, which is the norm of the difference vector.
diff_vector = p - q
length = linf_norm(diff_vector)

# Print the final calculation as a clear equation.
print("The length (L) of the segment is calculated as ||p - q||_inf.")
print(f"The vector for p is: {p}")
print(f"The vector for q is: {q}")
print(f"The difference vector p - q is: {diff_vector}\n")
print("The final equation is:")
# The following print statement shows all the numbers in the final equation.
print(f"L = ||({p[0]}, {p[1]}, {p[2]}) - ({q[0]}, {q[1]}, {q[2]})||_inf = ||({diff_vector[0]}, {diff_vector[1]}, {diff_vector[2]})||_inf = {length}")
