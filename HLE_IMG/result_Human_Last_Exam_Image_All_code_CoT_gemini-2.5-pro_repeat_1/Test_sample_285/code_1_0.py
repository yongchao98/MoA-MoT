import numpy as np

# This script analyzes triangulation B to check if it violates the Delaunay
# empty circle property.

# A violation occurs if the circumcircle of any triangle contains another data point.
# We suspect that in triangulation B, the circumcircle of the top-most triangle
# contains the top inner point.

# We assign plausible coordinates to the relevant points from image B.
# T = Top-most vertex
# TLo = Top-left outer vertex
# TRo = Top-right outer vertex
# Ti = Top inner vertex (the point we test)
T   = (0, 2.0)
TLo = (-2.0, 1.0)
TRo = (2.0, 1.0)
Ti  = (0, 0.5)

# The Delaunay in-circle test uses a determinant. For a triangle defined by
# points p1, p2, p3 in counter-clockwise (CCW) order, a fourth point p4
# is inside their circumcircle if the following determinant is positive.
#
# | p1.x  p1.y  p1.x^2+p1.y^2  1 |
# | p2.x  p2.y  p2.x^2+p2.y^2  1 |
# | p3.x  p3.y  p3.x^2+p3.y^2  1 |
# | p4.x  p4.y  p4.x^2+p4.y^2  1 |

# The order TLo -> TRo -> T is counter-clockwise.
p1 = TLo
p2 = TRo
p3 = T
p4 = Ti

# Construct the matrix for the determinant calculation.
# Each row corresponds to a point (x, y, x^2+y^2, 1).
matrix = np.array([
    [p1[0], p1[1], p1[0]**2 + p1[1]**2, 1],
    [p2[0], p2[1], p2[0]**2 + p2[1]**2, 1],
    [p3[0], p3[1], p3[0]**2 + p3[1]**2, 1],
    [p4[0], p4[1], p4[0]**2 + p4[1]**2, 1]
])

# Calculate the determinant of the matrix.
determinant = np.linalg.det(matrix)

print("Analysis of Triangulation B")
print("="*30)
print("We test if the top inner point Ti is inside the circumcircle of the triangle (TLo, TRo, T).")
print(f"Coordinates: TLo={p1}, TRo={p2}, T={p3}, Ti={p4}")
print("\nThe matrix for the in-circle test is:")
# The final equation is det(M) = result. We show the matrix M and the result.
for i in range(4):
    print(f"| {matrix[i, 0]:5.2f} {matrix[i, 1]:5.2f} {matrix[i, 2]:5.2f} {matrix[i, 3]:5.2f} |")

print(f"\nThe equation is det(M) = {determinant:.2f}")

print("\nConclusion:")
if determinant > 1e-9: # Use a small epsilon for floating point comparison
    print(f"The determinant is positive ({determinant:.2f}).")
    print("This confirms that the point Ti is INSIDE the circumcircle of the triangle (TLo, TRo, T).")
    print("Therefore, triangulation B violates the empty circle property.")
elif determinant < -1e-9:
    print(f"The determinant is negative ({determinant:.2f}).")
    print("The point Ti is OUTSIDE the circumcircle. This specific check does not show a violation.")
else:
    print(f"The determinant is zero.")
    print("The four points are co-circular, which is a degenerate case but doesn't violate the empty circle property.")
print("="*30)
<<<B>>>