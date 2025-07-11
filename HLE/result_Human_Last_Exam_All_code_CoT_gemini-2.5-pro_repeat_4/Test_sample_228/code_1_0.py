import numpy as np

def triangulate_point(p1, p2, P1, P2):
    """
    Triangulates a 3D point from two 2D observations and projection matrices
    using the Direct Linear Transform (DLT) method.
    """
    # Form the matrix A for AX=0
    A = np.array([
        p1[0] * P1[2, :] - P1[0, :],
        p1[1] * P1[2, :] - P1[1, :],
        p2[0] * P2[2, :] - P2[0, :],
        p2[1] * P2[2, :] - P2[1, :]
    ])

    # Solve AX=0 using SVD
    _, _, Vt = np.linalg.svd(A)
    X_homogeneous = Vt[-1]

    # Normalize to get Cartesian coordinates
    X_cartesian = X_homogeneous[:3] / X_homogeneous[3]
    return X_cartesian, X_homogeneous

# 1. Define the scene in a WORLD reference frame
# The true 3D point in the world frame
P_world = np.array([2.0, 3.0, 4.0])
P_world_h = np.append(P_world, 1) # Homogeneous coordinates

print(f"Original 3D point in WORLD frame: {P_world}")
print("-" * 40)

# 2. Define two cameras with their poses relative to the WORLD frame
# Camera intrinsics (same for both for simplicity)
K = np.array([[1000, 0, 320],
              [0, 1000, 240],
              [0, 0, 1]])

# Camera 1: At the origin of the world frame
R1 = np.identity(3)
t1 = np.zeros((3, 1))
T1_world_to_cam = np.hstack([R1, t1]) # Extrinsic matrix for Cam 1

# Camera 2: Shifted and rotated
R2 = np.array([[0.9, 0, 0.4],
               [0, 1, 0],
               [-0.4, 0, 0.9]]) # Simple rotation matrix
t2 = np.array([[-5], [0], [1]]) # Translation vector
T2_world_to_cam = np.hstack([R2, t2]) # Extrinsic matrix for Cam 2

# 3. Create full projection matrices (World -> Image)
P1 = K @ T1_world_to_cam
P2 = K @ T2_world_to_cam

# 4. Project the 3D world point onto each camera's 2D image plane
p1_h = P1 @ P_world_h
p1 = p1_h[:2] / p1_h[2]

p2_h = P2 @ P_world_h
p2 = p2_h[:2] / p2_h[2]

print(f"Projected 2D point in Camera 1: {p1}")
print(f"Projected 2D point in Camera 2: {p2}")
print("-" * 40)

# 5. Triangulate the 3D point using the 2D points and projection matrices.
# This process simulates triangulation using Plucker coordinates or other methods.
# The core idea is that it uses projection matrices defined from the world frame.
P_triangulated, P_triangulated_h = triangulate_point(p1, p2, P1, P2)

print("--- Triangulation Result ---")
print(f"Triangulated point coordinates: {P_triangulated}")
print(f"Note: The triangulated point is in the WORLD frame, same as the original point.")
print("-" * 40)

# 6. Show that a transformation is needed to get coordinates in a camera's frame.
# Let's find the coordinates of the triangulated point in Camera 2's reference frame.
P_in_cam2_frame_h = T2_world_to_cam @ P_triangulated_h
P_in_cam2_frame = P_in_cam2_frame_h[:3] / P_in_cam2_frame_h[3]

print("--- Frame Transformation ---")
print("To get the point's coordinates relative to a camera, we must apply a transformation.")
print(f"Triangulated point in CAMERA 2's frame: {P_in_cam2_frame.flatten()}")
print("This result is different from the world coordinates, proving that the triangulation")
print("does not directly yield a solution in the camera's frame.")
