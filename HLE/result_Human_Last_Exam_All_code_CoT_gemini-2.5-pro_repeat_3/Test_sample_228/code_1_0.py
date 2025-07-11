import numpy as np

def triangulate_plucker(l1, m1, l2, m2):
    """
    Finds the midpoint of the shortest segment between two skew lines
    represented by their Plücker coordinates in the same reference frame.
    L1 = (l1, m1), L2 = (l2, m2)
    """
    # Ensure direction vectors are normalized
    l1 = l1 / np.linalg.norm(l1)
    l2 = l2 / np.linalg.norm(l2)

    # Base points on each line (closest point to the origin)
    p1_base = np.cross(l1, m1)
    p2_base = np.cross(l2, m2)

    # Solve the 2x2 system to find the parameters for the closest points
    # A * [a, b]^T = B
    # where P1_closest = p1_base + a*l1 and P2_closest = p2_base + b*l2
    l1_dot_l2 = np.dot(l1, l2)
    A = np.array([[-1, l1_dot_l2],
                  [-l1_dot_l2, 1]])
    
    p_diff = p1_base - p2_base
    B = np.array([np.dot(p_diff, l1),
                  np.dot(p_diff, l2)])
    
    # Check if lines are parallel
    if np.abs(np.linalg.det(A)) < 1e-8:
        # For parallel lines, pick a simple average
        print("Warning: Lines are nearly parallel, triangulation may be unstable.")
        return (p1_base + p2_base) / 2

    # Solve for parameters a and b
    try:
        params = np.linalg.solve(A, B)
        a, b = params[0], params[1]
    except np.linalg.LinAlgError:
        print("Error: Could not solve for line parameters. Using pseudo-inverse.")
        params = np.linalg.pinv(A) @ B
        a, b = params[0], params[1]

    # Calculate the closest points on each line
    p1_closest = p1_base + a * l1
    p2_closest = p2_base + b * l2

    # The triangulated point is the midpoint of the shortest connecting segment
    triangulated_point = (p1_closest + p2_closest) / 2
    return triangulated_point

# --- Simulation Setup ---

# 1. Define a ground truth 3D point P in the world frame.
#    The world frame is assumed to be Camera 1's reference frame.
P_ground_truth = np.array([2.0, 3.0, 10.0])

# 2. Define Camera 1's pose (at the origin of the world frame).
#    The line of sight (ray) for Camera 1 goes from the origin towards P.
origin1 = np.array([0.0, 0.0, 0.0])
direction1 = P_ground_truth - origin1
direction1 = direction1 / np.linalg.norm(direction1) # Normalize

# 3. Define Camera 2's pose relative to Camera 1 (the world frame).
#    This transformation is REQUIRED to express Camera 2's ray in the world frame.
#    Let's say Camera 2 is shifted along the x-axis.
translation_cam2 = np.array([5.0, 0.0, 0.0])
origin2 = translation_cam2
# The line of sight for Camera 2 goes from its origin towards P.
direction2 = P_ground_truth - origin2
direction2 = direction2 / np.linalg.norm(direction2) # Normalize

# 4. Compute Plücker coordinates for both lines IN THE SAME (WORLD) FRAME.
#    L = (direction, moment = origin x direction)
l1 = direction1
m1 = np.cross(origin1, direction1) # This will be [0, 0, 0]

l2 = direction2
m2 = np.cross(origin2, direction2)

# 5. Triangulate the point using the Plücker coordinates.
#    Note: A small amount of noise could be added here to simulate a real scenario
#    where lines are skew, but for this perfect case, they will intersect.
P_triangulated = triangulate_plucker(l1, m1, l2, m2)

# 6. Print the result.
#    The result is in the world frame (Camera 1's frame).
#    This was only possible because we used the transformation (translation_cam2)
#    to define the second line in the same frame as the first.
print("This example demonstrates 3D point triangulation.")
print("The process requires expressing all lines of sight in a common reference frame.")
print("A transformation is needed to convert a ray from its native camera frame to this common frame.")
print("\n--- Results ---")
print(f"Original 3D Point:      ({P_ground_truth[0]:.4f}, {P_ground_truth[1]:.4f}, {P_ground_truth[2]:.4f})")

# Final output format as an equation
print("\nFinal Equation:")
print(f"Triangulated Point (X, Y, Z) = ({P_triangulated[0]:.4f}, {P_triangulated[1]:.4f}, {P_triangulated[2]:.4f})")