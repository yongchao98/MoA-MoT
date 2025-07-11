import numpy as np

def triangulate_skew_lines(p1, d1, p2, d2):
    """
    Finds the midpoint of the shortest segment between two skew lines.
    This midpoint is the optimal triangulated point.

    Args:
        p1 (np.ndarray): A point on the first line.
        d1 (np.ndarray): The direction vector of the first line.
        p2 (np.ndarray): A point on the second line.
        d2 (np.ndarray): The direction vector of the second line.

    Returns:
        np.ndarray: The triangulated 3D point.
        np.ndarray: The closest point on line 1.
        np.ndarray: The closest point on line 2.
    """
    # Normalize direction vectors for numerical stability
    d1 = d1 / np.linalg.norm(d1)
    d2 = d2 / np.linalg.norm(d2)

    # To find the closest points c1 and c2 on the lines, we solve a linear system.
    # The vector c1 - c2 must be perpendicular to both d1 and d2.
    # c1 = p1 + s * d1
    # c2 = p2 + t * d2
    # We solve for s and t using the conditions:
    # dot(c1 - c2, d1) = 0
    # dot(c1 - c2, d2) = 0

    # This results in a 2x2 linear system Ax = b for [s, -t]
    # [ dot(d1,d1)  dot(d1,d2) ] [  s ] = [ dot(p2-p1, d1) ]
    # [ dot(d2,d1)  dot(d2,d2) ] [ -t ] = [ dot(p2-p1, d2) ]

    d1d1 = np.dot(d1, d1)
    d2d2 = np.dot(d2, d2)
    d1d2 = np.dot(d1, d2)
    
    A = np.array([
        [d1d1, -d1d2],
        [d1d2, -d2d2]
    ])

    p_diff = p1 - p2
    b = np.array([-np.dot(p_diff, d1), -np.dot(p_diff, d2)])

    # Check if lines are parallel
    if np.abs(np.linalg.det(A)) < 1e-8:
        print("Lines are parallel or nearly parallel. Cannot compute a unique solution.")
        return None, None, None
        
    # Solve for s and t
    st = np.linalg.solve(A, b)
    s, t = st[0], st[1]

    # Calculate the closest points on each line
    c1 = p1 + s * d1
    c2 = p2 + t * d2

    # The triangulated point is the midpoint of the segment [c1, c2]
    triangulated_point = (c1 + c2) / 2.0
    
    return triangulated_point, c1, c2

# --- Main execution ---
# Define two skew lines in 3D space.
# These represent two back-projected rays from two cameras that fail to intersect.
# Line 1: Passes through (1, 0, 0), parallel to the Y-axis.
p1 = np.array([1.0, 0.0, 0.0])
d1 = np.array([0.0, 1.0, 0.0])

# Line 2: Passes through (0, 0, 1), parallel to the X-axis.
p2 = np.array([0.0, 0.0, 1.0])
d2 = np.array([1.0, 0.0, 0.0])

# Perform triangulation
P, c1, c2 = triangulate_skew_lines(p1, d1, p2, d2)

if P is not None:
    print("This script demonstrates triangulation for two non-intersecting (skew) lines.")
    print(f"Line 1 is defined by point p1={p1} and direction d1={d1}")
    print(f"Line 2 is defined by point p2={p2} and direction d2={d2}\n")
    print("Since the lines don't intersect, we find the points of closest approach on each line.")
    print(f"Closest point on Line 1 (c1): [{c1[0]:.2f}, {c1[1]:.2f}, {c1[2]:.2f}]")
    print(f"Closest point on Line 2 (c2): [{c2[0]:.2f}, {c2[1]:.2f}, {c2[2]:.2f}]\n")
    print("The final triangulated point (P) is the midpoint of the segment connecting c1 and c2.")
    
    # Print the final equation with all numbers
    print("Final Equation: P = (c1 + c2) / 2")
    c1_str = f"[{c1[0]:.2f}, {c1[1]:.2f}, {c1[2]:.2f}]"
    c2_str = f"[{c2[0]:.2f}, {c2[1]:.2f}, {c2[2]:.2f}]"
    p_str = f"[{P[0]:.2f}, {P[1]:.2f}, {P[2]:.2f}]"
    
    print(f"P = ({c1_str} + {c2_str}) / 2")
    print(f"Resulting Point P = {p_str}")
