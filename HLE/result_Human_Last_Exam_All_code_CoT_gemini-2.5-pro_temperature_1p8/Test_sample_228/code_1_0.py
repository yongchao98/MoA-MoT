import numpy as np

def triangulate_plucker(d1, m1, d2, m2):
    """
    Triangulates a 3D point from two skew 3D lines represented by
    Plucker coordinates.

    The method finds the midpoint of the shortest line segment connecting
    the two skew lines. The lines and the resulting point are all in the
    same coordinate frame.

    Args:
        d1, m1: Direction and moment vectors for the first line (L1).
        d2, m2: Direction and moment vectors for the second line (L2).

    Returns:
        The triangulated 3D point, and the two closest points on each line.
    """
    d1, m1 = np.asarray(d1), np.asarray(m1)
    d2, m2 = np.asarray(d2), np.asarray(m2)

    # Check for parallel lines, which is a degenerate case.
    if np.allclose(np.cross(d1, d2), 0):
        print("Warning: Lines are parallel. Solution is not uniquely defined.")
        return None, None, None

    # A point on a line is given by A = (d x m) / ||d||^2 + k*d.
    # We find the base points A1, A2 (closest points on each line to the origin).
    A1 = np.cross(d1, m1) / np.dot(d1, d1)
    A2 = np.cross(d2, m2) / np.dot(d2, d2)

    # We need to find parameters t and s for points P1 = A1 + t*d1 and
    # P2 = A2 + s*d2 such that the vector P1-P2 is perpendicular to both d1 and d2.
    # This leads to a 2x2 linear system for t and s:
    # (A1 + t*d1 - (A2 + s*d2)) · d1 = 0
    # (A1 + t*d1 - (A2 + s*d2)) · d2 = 0
    #
    # Rearranging gives:
    # t*(d1·d1) - s*(d1·d2) = (A2-A1)·d1
    # t*(d1·d2) - s*(d2·d2) = (A2-A1)·d2

    # Setup the 2x2 matrix M and vector b
    M = np.array([
        [np.dot(d1, d1), -np.dot(d1, d2)],
        [np.dot(d1, d2), -np.dot(d2, d2)]
    ])
    b = np.array([np.dot(A2 - A1, d1), np.dot(A2 - A1, d2)])

    # Solve the system M*x = b for x = [t, s]
    try:
        params = np.linalg.solve(M, b)
        t, s = params[0], params[1]
    except np.linalg.LinAlgError:
        print("Error: Could not solve for line parameters.")
        return None, None, None

    # Calculate the closest points on each line
    P1 = A1 + t * d1
    P2 = A2 + s * d2

    # The triangulated point is the midpoint of the segment connecting P1 and P2
    P_triangulated = (P1 + P2) / 2.0

    return P_triangulated, P1, P2

# --- Main execution ---
# Define two skew lines in a common reference frame using points on the lines.
# Plücker coordinates are L = (d, m) where d is direction and m = p1 x p2.

# Line 1 (L1): A line parallel to the Y-axis at x=1, z=0.
P_a1 = np.array([1.0, 0.0, 0.0])
P_b1 = np.array([1.0, 1.0, 0.0])
d1 = P_b1 - P_a1
m1 = np.cross(P_a1, P_b1)

# Line 2 (L2): A line parallel to the X-axis at y=0, z=1.
P_a2 = np.array([0.0, 0.0, 1.0])
P_b2 = np.array([1.0, 0.0, 1.0])
d2 = P_b2 - P_a2
m2 = np.cross(P_a2, P_b2)

# Check if lines intersect using the reciprocal product (d1·m2 + d2·m1)
reciprocal_product = np.dot(d1, m2) + np.dot(d2, m1)

print("--- Line Definitions (in a common reference frame) ---")
print(f"L1 Direction (d1): {d1}")
print(f"L1 Moment (m1):    {m1}")
print(f"L2 Direction (d2): {d2}")
print(f"L2 Moment (m2):    {m2}")
print(f"\nReciprocal product (d1·m2 + d2·m1): {reciprocal_product:.4f}")

if np.isclose(reciprocal_product, 0):
    print("The lines intersect geometrically.")
else:
    print("The lines are skew (do not intersect). This is the key limitation.")

print("\n--- Triangulation of Skew Lines ---")
P, P1_closest, P2_closest = triangulate_plucker(d1, m1, d2, m2)

if P is not None:
    print(f"Point on L1 closest to L2 (P1):   {np.round(P1_closest, 3)}")
    print(f"Point on L2 closest to L1 (P2):   {np.round(P2_closest, 3)}")
    print(f"Shortest distance between lines:  {np.linalg.norm(P1_closest - P2_closest):.4f}")
    
    print("\nFinal Triangulated Point (midpoint of shortest segment):")
    # We print the final equation showing all the numbers involved.
    print(f"P = ( P1 + P2 ) / 2.0")
    print(f"P = ( [{P1_closest[0]:.3f}, {P1_closest[1]:.3f}, {P1_closest[2]:.3f}] + [{P2_closest[0]:.3f}, {P2_closest[1]:.3f}, {P2_closest[2]:.3f}] ) / 2.0")
    print(f"P = [ {P[0]:.3f}, {P[1]:.3f}, {P[2]:.3f} ]")
