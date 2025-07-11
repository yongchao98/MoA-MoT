import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    This matrix represents the cross-product operation, i.e., skew(v) @ u = v x u.
    """
    if len(v) != 3:
        raise ValueError("Input vector must have 3 elements")
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

# Plan:
# 1. Define a 3D point in the world frame.
# 2. Define two cameras. Cam1 is at the world origin. Cam2 is translated.
# 3. For each camera, determine the line of sight to the 3D point.
# 4. Represent these lines using Plücker coordinates in a COMMON reference frame (the world frame).
#    This step demonstrates the need for transformation for Cam2's line.
# 5. Triangulate the 3D point using the two Plücker lines in the common frame.
# 6. Compare the result with the original point.

# 1. Define a 3D point in the world coordinate system
P_world = np.array([2, 3, 10])

# 2. Define camera poses (Rotation | translation) relative to the world frame
# Camera 1 is at the origin of the world frame
R1 = np.identity(3)
t1 = np.array([0, 0, 0])
C1 = t1 # Camera center 1

# Camera 2 is translated along the x-axis
R2 = np.identity(3)
t2 = np.array([5, 0, 0])
C2 = t2 # Camera center 2

print("This script demonstrates that a transformation is necessary for triangulation with Plücker coordinates.")
print("We will triangulate a point from two camera views.\n")

# 3. Determine the line of sight for each camera IN THE WORLD FRAME
# Line of sight is a 3D line passing through the camera center and the 3D point.

# For Line 1 (from Camera 1)
# Direction vector in world frame
d1_world = (P_world - C1)
d1_world = d1_world / np.linalg.norm(d1_world)
# Moment vector in world frame (m = C x d)
m1_world = np.cross(C1, d1_world)

# For Line 2 (from Camera 2)
# To get the Plücker coordinates in the world frame, we must use world coordinates.
# We CANNOT use the line's local camera coordinates directly.
# This transformation from a local camera frame to the world frame is the key step.
# Here, we compute it directly in the world frame for simplicity.
d2_world = (P_world - C2)
d2_world = d2_world / np.linalg.norm(d2_world)
# Moment vector in world frame (m = C x d)
m2_world = np.cross(C2, d2_world)

print("--- Final Equation Components ---")
print("The triangulation finds a point P that lies on two lines, L1 and L2.")
print("This can be expressed as a system of linear equations derived from m = P x d for each line.")
print("L1 direction (d1):", d1_world)
print("L1 moment (m1):  ", m1_world)
print("L2 direction (d2):", d2_world)
print("L2 moment (m2):  ", m2_world)
print("\nNote: The coordinates for L2 are computed relative to the world frame, not L2's local frame.")
print("This implicitly performs the necessary transformation.\n")


# 5. Triangulate the 3D point using the two Plücker lines
# We need to solve for a point P such that:
# m1 = P x d1  =>  [d1]_x * P = m1
# m2 = P x d2  =>  [d2]_x * P = m2
# where [d]_x is the skew-symmetric matrix of d.
# We form an overdetermined system Ax = b and solve using the pseudo-inverse.

A = np.vstack([skew(d1_world), skew(d2_world)])
b = np.hstack([m1_world, m2_world])

# Solve for P using the pseudo-inverse: P = (A^T * A)^-1 * A^T * b
P_triangulated = np.linalg.pinv(A) @ b

# 6. Compare results
print("--- Triangulation Result ---")
print("Original 3D point P:    ", P_world)
print("Triangulated 3D point P:", P_triangulated.round(decimals=5))

<<<E>>>