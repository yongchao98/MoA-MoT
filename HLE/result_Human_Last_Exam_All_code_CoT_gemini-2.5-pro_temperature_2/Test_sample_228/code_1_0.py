import numpy as np

def triangulate_point(p1, p2, P1, P2):
    """
    Triangulates a 3D point from two 2D observations and camera matrices
    using the Direct Linear Transform (DLT) method.
    """
    A = np.zeros((4, 4))
    A[0] = p1[0] * P1[2, :] - P1[0, :]
    A[1] = p1[1] * P1[2, :] - P1[1, :]
    A[2] = p2[0] * P2[2, :] - P2[0, :]
    A[3] = p2[1] * P2[2, :] - P2[1, :]
    
    # Solve the system AX=0 using Singular Value Decomposition (SVD)
    _, _, Vh = np.linalg.svd(A)
    # The solution is the last column of V, which is the last row of Vh
    X_homogeneous = Vh[-1]
    
    # Dehomogenize to get 3D coordinates
    X = X_homogeneous[:3] / X_homogeneous[3]
    return X

# 1. Setup the 3D Scene
# Define a true 3D point in the WORLD coordinate system
P_world = np.array([2, 3, 10])

# Define a camera intrinsic matrix (same for both cameras for simplicity)
K = np.array([[800, 0,   320],
              [0,   800, 240],
              [0,   0,   1]])

# Camera 1: At the origin of the world coordinate system
R1 = np.eye(3)
t1 = np.zeros((3, 1))
P1 = K @ np.hstack((R1, t1)) # Projection matrix for Camera 1

# Camera 2: Rotated and translated in the world
# Rotation of -15 degrees around the y-axis
theta = np.deg2rad(-15)
R2 = np.array([[np.cos(theta),  0, np.sin(theta)],
               [0,              1, 0],
               [-np.sin(theta), 0, np.cos(theta)]])
t2 = np.array([[-5], [0], [0]]) # Translation vector for Camera 2
P2 = K @ np.hstack((R2, t2)) # Projection matrix for Camera 2

# 2. Simulate Image Formation
# Project the 3D point onto the image planes
p1_homogeneous = P1 @ np.append(P_world, 1)
p1 = p1_homogeneous[:2] / p1_homogeneous[2]

p2_homogeneous = P2 @ np.append(P_world, 1)
p2 = p2_homogeneous[:2] / p2_homogeneous[2]

print("--- Setup ---")
print(f"Original 3D Point in World Frame (X_world):")
print(f"X={P_world[0]}, Y={P_world[1]}, Z={P_world[2]}\n")

# 3. Perform Triangulation
# This process uses p1, p2, P1, and P2, where P1 and P2 map from the WORLD frame.
# Therefore, the output must also be in the WORLD frame.
# This logic is the same for any triangulation method, including those using Plücker coordinates.
P_triangulated_world = triangulate_point(p1, p2, P1, P2)

print("--- Triangulation Result ---")
print("Triangulation produces the point in the common 'World' reference frame.")
print("Triangulated Point in World Frame (X_triangulated_world):")
# Output each number in the final equation/result
print(f"X={P_triangulated_world[0]:.4f}, Y={P_triangulated_world[1]:.4f}, Z={P_triangulated_world[2]:.4f}\n")


# 4. Analyze the Result and Transform
# To get the point in a specific CAMERA's reference frame, a transformation is needed.
# Let's find the coordinates of the point in Camera 2's frame.
# We apply the extrinsic transformation of Camera 2: P_cam = R * P_world + t
P_in_cam2_frame = (R2 @ P_triangulated_world).reshape(3, 1) + t2

print("--- Transformation to Camera Frame ---")
print("To get the coordinates in a specific camera's frame, we must apply a transformation.")
print(f"The equation is: X_camera2 = R2 * X_world + t2")
print("\nPoint Coordinates in Camera 2's Reference Frame:")
# Output each number in the final equation/result
print(f"X={P_in_cam2_frame[0][0]:.4f}, Y={P_in_cam2_frame[1][0]:.4f}, Z={P_in_cam2_frame[2][0]:.4f}")
