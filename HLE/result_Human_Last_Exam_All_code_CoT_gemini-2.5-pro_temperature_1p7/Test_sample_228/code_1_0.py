import numpy as np

def triangulate_point(p1, p2, K1, M1, K2, M2):
    """
    Triangulates a 3D point from two 2D correspondences and camera matrices.
    This implementation solves the linear system for the point of closest approach.

    Args:
        p1 (np.array): 2D point in image 1 (shape 2,).
        p2 (np.array): 2D point in image 2 (shape 2,).
        K1 (np.array): Camera 1 intrinsic matrix (3x3).
        M1 (np.array): Camera 1 extrinsic matrix [R|t] (3x4).
        K2 (np.array): Camera 2 intrinsic matrix (3x3).
        M2 (np.array): Camera 2 extrinsic matrix [R|t] (3x4).

    Returns:
        np.array: The triangulated 3D point in the world frame.
    """
    # 1. Define lines in the world frame
    # Camera centers are the origins of the rays.
    R1 = M1[:, :3]
    t1 = M1[:, 3]
    # C = -R' * t
    C1 = -R1.T @ t1

    R2 = M2[:, :3]
    t2 = M2[:, 3]
    C2 = -R2.T @ t2

    # Directions of the rays.
    # Back-project 2D point to a 3D ray in the camera's frame, then rotate to world frame.
    p1_h = np.array([p1[0], p1[1], 1.0])
    d1_cam = np.linalg.inv(K1) @ p1_h
    d1_world = R1.T @ d1_cam
    d1_world /= np.linalg.norm(d1_world)

    print("Step 1: The line from Camera 1 is defined by a point (its center) and a direction vector, both in the world frame.")
    print(f"Camera 1 Center (C1): [{C1[0]:.2f}, {C1[1]:.2f}, {C1[2]:.2f}]")
    print(f"Ray 1 Direction (d1_world): [{d1_world[0]:.2f}, {d1_world[1]:.2f}, {d1_world[2]:.2f}]\n")


    # For the second camera, we MUST transform its local ray into the world frame.
    # This is the key transformation the problem asks about.
    p2_h = np.array([p2[0], p2[1], 1.0])
    d2_cam = np.linalg.inv(K2) @ p2_h
    # Transform direction vector from camera 2's frame to the world frame using its rotation matrix.
    d2_world = R2.T @ d2_cam
    d2_world /= np.linalg.norm(d2_world)
    print("Step 2: A transformation is needed for Camera 2.")
    print("The ray from Camera 2 must be transformed from its local camera frame into the common world frame.")
    print(f"Camera 2 Center (C2): [{C2[0]:.2f}, {C2[1]:.2f}, {C2[2]:.2f}]")
    print(f"Ray 2 Direction (d2_world): [{d2_world[0]:.2f}, {d2_world[1]:.2f}, {d2_world[2]:.2f}]\n")


    # 2. Find the point of closest approach between the two lines in the world frame.
    # We solve the system `(C1 + s1*d1) - (C2 + s2*d2)` is perpendicular to d1 and d2.
    # This gives a 2x2 linear system for s1 and s2.
    A = np.array([
        [np.dot(d1_world, d1_world), -np.dot(d1_world, d2_world)],
        [np.dot(d1_world, d2_world), -np.dot(d2_world, d2_world)]
    ])
    b = np.array([
        np.dot(C2 - C1, d1_world),
        np.dot(C2 - C1, d2_world)
    ])
    
    # Check if lines are parallel
    if np.abs(np.linalg.det(A)) < 1e-8:
        print("Warning: Lines are nearly parallel. Triangulation may be unstable.")
        # Simplified handling for parallel lines: return midpoint of camera centers projection
        midpoint_centers = (C1 + C2) / 2
        return midpoint_centers
        
    s, _ = np.linalg.solve(A, b)
    s1 = s[0]
    s2 = s[1]

    # The two points on each line that are closest to each other
    P1 = C1 + s1 * d1_world
    P2 = C2 + s2 * d2_world

    # The triangulated point is the midpoint of the segment connecting P1 and P2
    P_world_recon = (P1 + P2) / 2
    return P_world_recon


# --- Simulation Setup ---
# Intrinsics (assuming same for both cameras)
K = np.array([
    [800, 0, 320],
    [0, 800, 240],
    [0, 0, 1]
])

# Extrinsics
# Camera 1 is at the world origin
R1 = np.eye(3)
t1 = np.zeros(3)
M1 = np.c_[R1, t1]

# Camera 2 is translated along the x-axis and slightly rotated
# Rotation matrix for a 15-degree rotation around y-axis
angle = np.deg2rad(-15)
R2 = np.array([
    [np.cos(angle), 0, np.sin(angle)],
    [0, 1, 0],
    [-np.sin(angle), 0, np.cos(angle)]
])
t2 = np.array([-2, 0.1, 0.5])
M2 = np.c_[R2, t2]

# Define a true 3D point in the world
P_world_true = np.array([0.5, 0.2, 5.0])

# Project the 3D point into both cameras to get 2D image points
P_world_h = np.append(P_world_true, 1) # Homogeneous coordinates
p1_h = K @ M1 @ P_world_h
p1 = p1_h[:2] / p1_h[2]

p2_h = K @ M2 @ P_world_h
p2 = p2_h[:2] / p2_h[2]

# --- Run Triangulation ---
print("This script demonstrates that a transformation is necessary for triangulation.")
print("Even if using Plucker coordinates, lines from different cameras must first be put in a common reference frame.\n")
P_world_reconstructed = triangulate_point(p1, p2, K, M1, K, M2)

print("Step 3: With lines in a common frame, we can find the 3D point.")
print("The final triangulated point is calculated as the midpoint of the shortest segment connecting the two rays.\n")
print("Original 3D point:")
print(f"[{P_world_true[0]:.4f}, {P_world_true[1]:.4f}, {P_world_true[2]:.4f}]")
print("\nReconstructed 3D point after transformation and triangulation:")
# Final output formatting: X, Y, Z
print(f"[{P_world_reconstructed[0]:.4f}, {P_world_reconstructed[1]:.4f}, {P_world_reconstructed[2]:.4f}]")