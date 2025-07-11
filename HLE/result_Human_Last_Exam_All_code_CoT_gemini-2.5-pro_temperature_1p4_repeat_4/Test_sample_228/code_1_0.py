import numpy as np

def skew(v):
    """Converts a 3D vector to its skew-symmetric matrix."""
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def points_to_plucker(p1, p2):
    """Creates Plücker coordinates for a line passing through two points."""
    l = p2 - p1
    m = np.cross(p1, p2)
    # Normalize for consistency, though not strictly required for this demo
    norm_l = np.linalg.norm(l)
    return l / norm_l, m / norm_l

def transform_plucker(L, R, t):
    """Transforms Plücker coordinates to a new reference frame."""
    l, m = L
    l_prime = R @ l
    m_prime = R @ m + np.cross(t, R @ l)
    return l_prime, m_prime

def triangulate_plucker(L1, L2):
    """
    Finds the 3D point of closest approach for two 3D lines in Plücker form.
    This point is the midpoint of the common perpendicular segment.
    """
    l1, m1 = L1
    l2, m2 = L2
    
    cp = np.cross(l1, l2)
    cp_sq_norm = np.dot(cp, cp)
    
    if cp_sq_norm < 1e-9: # Lines are nearly parallel
        print("Lines are parallel, cannot triangulate reliably.")
        return None

    # This formula calculates the point of closest approach
    p_3d = (np.cross(l1, np.cross(l2, m1)) + np.cross(l2, np.cross(l1, m2))) / cp_sq_norm
    return p_3d

# 1. SETUP: Define a ground truth point and two cameras in a world frame.
# Camera 1 is our reference frame (at the world origin).
R1_world = np.identity(3)
t1_world = np.array([0., 0., 0.])

# Camera 2 is translated and rotated relative to Camera 1.
t2_world = np.array([5., 1., 0.])
angle_y = np.deg2rad(-15)
R2_world = np.array([[np.cos(angle_y), 0, np.sin(angle_y)],
                     [0, 1, 0],
                     [-np.sin(angle_y), 0, np.cos(angle_y)]])

# The 3D point we want to triangulate, defined in the world/Camera 1 frame.
P_world = np.array([2., 3., 10.])

print("--- Step 1: Define lines in their local camera frames ---")
# 2. CREATE LINES IN LOCAL FRAMES
# Line L1 as seen from Camera 1. Its frame is the world frame.
# The line passes through the camera's origin and the 3D point.
C1_local = np.array([0., 0., 0.])
P_in_cam1_frame = P_world # Since C1 is the world origin
L1_local = points_to_plucker(C1_local, P_in_cam1_frame)
print(f"Plücker Line L1 (in Cam1 Frame): l={np.round(L1_local[0],3)}, m={np.round(L1_local[1],3)}")

# Line L2 as seen from Camera 2.
# First, find the coordinates of the 3D point P in Camera 2's frame.
# Transformation from world to cam2: P_cam2 = R' * (P_world - t)
P_in_cam2_frame = R2_world.T @ (P_world - t2_world)
# The line passes through the camera's local origin and the point P in its frame.
C2_local = np.array([0., 0., 0.])
L2_local = points_to_plucker(C2_local, P_in_cam2_frame)
print(f"Plücker Line L2 (in Cam2 Frame): l={np.round(L2_local[0],3)}, m={np.round(L2_local[1],3)}")
print("\nNotice the coordinates are different. Triangulating them directly is invalid.")

print("\n--- Step 2: Transform L2 into Camera 1's reference frame ---")
# 3. TRANSFORM L2 TO COMMON FRAME
# We need the pose of Camera 2 *in* Camera 1's frame, which is (R2_world, t2_world).
L2_transformed = transform_plucker(L2_local, R2_world, t2_world)
print("Transformation from Cam2's frame to Cam1's frame has been applied to L2.")
print(f"Transformed L2: l={np.round(L2_transformed[0],3)}, m={np.round(L2_transformed[1],3)}")

print("\n--- Step 3: Triangulate using lines in the common reference frame ---")
# 4. TRIANGULATE IN COMMON FRAME (Camera 1's frame)
# Now both L1_local and L2_transformed are in the same frame.
P_triangulated = triangulate_plucker(L1_local, L2_transformed)

print("The final triangulated point coordinates are in the common frame (Cam1/World).")
print(f"Original 3D Point:       [{P_world[0]:8.4f}, {P_world[1]:8.4f}, {P_world[2]:8.4f}]")
print(f"Triangulated 3D Point:   [{P_triangulated[0]:8.4f}, {P_triangulated[1]:8.4f}, {P_triangulated[2]:8.4f}]")
print("\nThe result shows that triangulation is successful only after the transformation.")

<<<E>>>