import numpy as np

def triangulate_point(P1, P2, p1, p2):
    """
    Triangulates a 3D point from two camera views using the Direct Linear Transform (DLT).

    Args:
        P1 (3x4 array): Projection matrix for camera 1.
        P2 (3x4 array): Projection matrix for camera 2.
        p1 (2x1 array): Homogeneous coordinates of the point in image 1.
        p2 (2x1 array): Homogeneous coordinates of the point in image 2.

    Returns:
        X (4x1 array): The triangulated 3D point in homogeneous coordinates (in the world frame).
    """
    # Form the matrix A for the linear system Ax = 0
    A = np.array([
        p1[0] * P1[2, :] - P1[0, :],
        p1[1] * P1[2, :] - P1[1, :],
        p2[0] * P2[2, :] - P2[0, :],
        p2[1] * P2[2, :] - P2[1, :]
    ])

    # Solve Ax = 0 using Singular Value Decomposition (SVD)
    _, _, Vh = np.linalg.svd(A)
    X = Vh[-1, :]
    
    # Normalize to make the last component 1
    return X / X[3]

# 1. Define a ground truth 3D point in the WORLD reference frame.
P_world = np.array([2.0, 3.0, 4.0, 1.0])

# 2. Define two cameras
# Camera 1 is at the origin of the world frame
R1 = np.eye(3)
t1 = np.zeros((3, 1))
P1 = np.hstack((R1, t1)) # Projection matrix for Camera 1

# Camera 2 is translated and rotated
# Rotation around Y axis by -10 degrees
angle = np.deg2rad(-10)
R2 = np.array([
    [np.cos(angle), 0, np.sin(angle)],
    [0, 1, 0],
    [-np.sin(angle), 0, np.cos(angle)]
])
# Translation
t2 = np.array([[-5.0], [0.0], [1.0]])
P2 = np.hstack((R2, t2)) # Projection matrix for Camera 2

print("--- Setup ---")
print(f"Original 3D point in WORLD frame: {P_world[:3]}")
print(f"Camera 2 is positioned relative to the world frame.")
print("-" * 20 + "\n")

# 3. Project the 3D point onto each camera's image plane to get 2D points
# (This simulates the process of observing the point)
p1_proj = P1 @ P_world
p1 = p1_proj[:2] / p1_proj[2]

p2_proj = P2 @ P_world
p2 = p2_proj[:2] / p2_proj[2]

print("--- Triangulation ---")
print("Using the 2D projections from each camera to triangulate the 3D point.")

# 4. Triangulate the 3D point using the projection matrices and 2D points
P_triangulated_world = triangulate_point(P1, P2, p1, p2)

print(f"Triangulated point coordinates: {P_triangulated_world[:3]}")
print("Observation: The result of triangulation is a point in the WORLD reference frame.")
print("It matches our original world point, not a point relative to either camera.\n")

# 5. Explicitly transform the world point into Camera 2's reference frame
# This shows the necessary extra step.
P_world_3d = P_world[:3].reshape(3, 1)
P_cam2 = R2 @ P_world_3d + t2

print("--- Transformation to Camera Frame ---")
print("To get the point's coordinates in Camera 2's frame, we must apply a transformation:")
print("P_cam2 = R2 * P_world + t2\n")
print(f"Original point coordinates in CAMERA 2 frame: \n{P_cam2.flatten()}\n")

# We can do the same for the triangulated point to verify
P_triangulated_3d = P_triangulated_world[:3].reshape(3, 1)
P_triangulated_cam2 = R2 @ P_triangulated_3d + t2
print(f"Triangulated point coordinates in CAMERA 2 frame: \n{P_triangulated_cam2.flatten()}\n")
print("Conclusion: The triangulation yields a world-frame point. A transformation is needed to get the camera-frame point.")
