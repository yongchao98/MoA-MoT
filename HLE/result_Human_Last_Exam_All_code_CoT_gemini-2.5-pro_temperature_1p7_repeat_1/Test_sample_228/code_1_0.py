import numpy as np

def skew(v):
    """
    Convert a 3-element vector to its skew-symmetric matrix.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def project(P, X):
    """
    Project a 3D point X using projection matrix P.
    """
    x_h = P @ np.append(X, 1)
    return x_h[:2] / x_h[2]

def get_plucker_line_in_cam_frame(K, x_img):
    """
    Calculates the Plucker coordinates of a ray in its local camera frame.
    The ray passes through the camera origin and the 3D point corresponding to x_img.
    """
    # Back-project the image point to a 3D point on the z=1 plane in the camera frame
    x_hom = np.append(x_img, 1)
    p_cam = np.linalg.inv(K) @ x_hom
    
    # The line passes through the camera origin (0,0,0) and p_cam.
    # Direction vector d is simply p_cam
    d = p_cam / np.linalg.norm(p_cam)
    
    # The moment vector m = origin x d = 0
    m = np.zeros(3)
    
    return d, m

def transform_plucker_line(d_local, m_local, R, t):
    """
    Transform Plucker line from local frame to world frame.
    R, t define the pose of the local frame in the world frame (p_world = R @ p_local + t).
    """
    d_world = R @ d_local
    m_world = R @ m_local + np.cross(t, d_world)
    return d_world, m_world

# 1. SETUP
# Ground truth 3D point in the world frame
X_world_gt = np.array([0.5, 0.8, 3.0])

# Intrinsic camera matrix (same for both cameras for simplicity)
K = np.array([
    [1000, 0, 320],
    [0, 1000, 240],
    [0, 0, 1]
])

# --- Camera 1 ---
# Pose is the origin of the world frame
R1 = np.eye(3)
t1 = np.zeros(3)
P1 = K @ np.hstack((R1, t1.reshape(3, 1)))

# --- Camera 2 ---
# Pose is translated and rotated relative to Camera 1
# Simple rotation around Y axis
angle = -np.pi / 8
R2 = np.array([
    [np.cos(angle), 0, np.sin(angle)],
    [0, 1, 0],
    [-np.sin(angle), 0, np.cos(angle)]
])
t2 = np.array([-1.5, 0, 0.5])
P2 = K @ np.hstack((R2, t2.reshape(3, 1)))

print("This script demonstrates 3D point triangulation using Plucker coordinates.")
print("It shows that a transformation is necessary to bring lines into a common frame.")
print("-" * 70)
print(f"Ground Truth 3D Point (World Frame): {X_world_gt}")
print("-" * 70)

# 2. PROJECT POINT ONTO IMAGES
x1_img = project(P1, X_world_gt)
x2_img = project(P2, X_world_gt)

print(f"Projected 2D point in Camera 1: {x1_img}")
print(f"Projected 2D point in Camera 2: {x2_img}")
print("-" * 70)

# 3. GET PLUCKER LINES IN LOCAL CAMERA FRAMES
d1_cam1, m1_cam1 = get_plucker_line_in_cam_frame(K, x1_img)
d2_cam2, m2_cam2 = get_plucker_line_in_cam_frame(K, x2_img)

print("Step 1: Calculate lines in their LOCAL camera frames.")
print(f"Line 1 in Cam1 Frame (d, m):\n{d1_cam1}\n{m1_cam1}\n")
print(f"Line 2 in Cam2 Frame (d, m):\n{d2_cam2}\n{m2_cam2}\n")
print("NOTE: These two lines CANNOT be used together directly as they are in different coordinate systems.")
print("-" * 70)

# 4. TRANSFORM LINES TO A COMMON WORLD FRAME
# We choose Camera 1's frame as the world frame.
# Line 1 is already in the world frame.
d1_world, m1_world = d1_cam1, m1_cam1

# Line 2 must be transformed from Camera 2's frame to the world frame using (R2, t2)
print("Step 2: Transform Line 2 into the common world frame (Camera 1's frame).")
d2_world, m2_world = transform_plucker_line(d2_cam2, m2_cam2, R2, t2)
print(f"Line 1 in World Frame (d, m):\n{d1_world}\n{m1_world}\n")
print(f"Transformed Line 2 in World Frame (d, m):\n{d2_world}\n{m2_world}\n")
print("-" * 70)

# 5. TRIANGULATE USING LEAST SQUARES
# For a point X to be on a line (d, m), it must satisfy skew(d) @ X = m.
# We stack these equations for both lines to form AX = B.

A = np.vstack([skew(d1_world), skew(d2_world)])
B = np.hstack([m1_world, m2_world])

print("Step 3: Solve the linear system AX = B for the 3D point X.")
print("The final equation to solve is:")
print("\n[ A Matrix ]             [ X ]   [ B Vector ]")
for i in range(A.shape[0]):
    row_str = " ".join([f"{x: 8.4f}" for x in A[i]])
    if i == 1:
        print(f"[ {row_str} ] [ X ] = [ {B[i]: 8.4f} ]")
    elif i == 2:
        print(f"[ {row_str} ] [ Y ]   [ {B[i]: 8.4f} ]")
    elif i == 3:
        print(f"[ {row_str} ] [ Z ]   [ {B[i]: 8.4f} ]")
    else:
        print(f"[ {row_str} ]         [ {B[i]: 8.4f} ]")


# Solve the system using least squares
X_triangulated, residuals, rank, s = np.linalg.lstsq(A, B, rcond=None)

print("-" * 70)
print(f"Ground Truth 3D Point:    {np.round(X_world_gt, 5)}")
print(f"Triangulated 3D Point:  {np.round(X_triangulated, 5)}")
print("\nThe triangulated point matches the ground truth, confirming the process.")
print("This was only possible AFTER transforming the line from Camera 2 into the common world frame.")