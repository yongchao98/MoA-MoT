import numpy as np

def triangulate_point(p1, p2, K, pose1, pose2):
    """
    Triangulates a 3D point from two 2D observations and camera parameters.
    This function demonstrates that line representations from each camera must be
    transformed into a common coordinate frame (the world frame) before triangulation.
    """
    # Decompose camera poses into rotation (R) and translation (t)
    R1, t1 = pose1[:, :3], pose1[:, 3]
    R2, t2 = pose2[:, :3], pose2[:, 3]

    # --- Step 1: Back-project 2D points to rays in local camera coordinates ---
    # The ray is defined by its origin (0,0,0) and a direction vector.
    # We find the direction by converting the 2D image point to a 3D point on the
    # image plane (at Z=1) and then un-projecting it using the inverse intrinsic matrix.
    p1_homogeneous = np.array([p1[0], p1[1], 1.0])
    p2_homogeneous = np.array([p2[0], p2[1], 1.0])

    K_inv = np.linalg.inv(K)
    d_cam1 = K_inv @ p1_homogeneous
    d_cam1 = d_cam1 / np.linalg.norm(d_cam1) # Normalize direction vector

    d_cam2 = K_inv @ p2_homogeneous
    d_cam2 = d_cam2 / np.linalg.norm(d_cam2) # Normalize direction vector

    print("--- Ray representation in LOCAL camera frames ---")
    print(f"Ray 1 direction (in Cam1 frame): {np.round(d_cam1, 3)}")
    print(f"Ray 2 direction (in Cam2 frame): {np.round(d_cam2, 3)}")
    print("\nThese rays cannot be directly used together as they are in different coordinate systems.")

    # --- Step 2: Transform rays into the common WORLD coordinate frame ---
    # This is the essential transformation that the problem asks about.
    # The origin of each ray is the camera's center in the world.
    # The direction of each ray is its camera-frame direction rotated into the world frame.
    o1_world = t1
    d1_world = R1 @ d_cam1

    o2_world = t2
    d2_world = R2 @ d_cam2

    print("\n--- Ray representation in the common WORLD frame ---")
    print(f"Ray 1 origin (world): {np.round(o1_world, 3)}")
    print(f"Ray 1 direction (world): {np.round(d1_world, 3)}")
    print(f"Ray 2 origin (world): {np.round(o2_world, 3)}")
    print(f"Ray 2 direction (world): {np.round(d2_world, 3)}")

    # --- Step 3: Triangulate by finding the closest point between the two rays ---
    # We solve a linear system to find the parameters 'a' and 'b' for the points
    # P1 = o1 + a*d1 and P2 = o2 + b*d2 on each line that are closest to each other.
    # The system is A * x = B, where x = [a, b]^T
    d1_dot_d1 = np.dot(d1_world, d1_world)
    d1_dot_d2 = np.dot(d1_world, d2_world)
    d2_dot_d2 = np.dot(d2_world, d2_world)
    o_diff = o2_world - o1_world

    A = np.array([
        [d1_dot_d1, -d1_dot_d2],
        [d1_dot_d2, -d2_dot_d2]
    ])

    B = np.array([
        np.dot(d1_world, o_diff),
        np.dot(d2_world, o_diff)
    ])

    print("\n--- Solving Linear System Ax = B for Triangulation ---")
    print(f"Matrix A = [[{A[0,0]:.3f}, {A[0,1]:.3f}], [{A[1,0]:.3f}, {A[1,1]:.3f}]]")
    print(f"Vector B = [{B[0]:.3f}, {B[1]:.3f}]")

    # Solve for scalars a and b
    try:
        a, b = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("Lines are parallel, cannot triangulate.")
        return None

    # The triangulated point is the midpoint of the shortest segment connecting the lines
    P1 = o1_world + a * d1_world
    P2 = o2_world - b * d2_world # Note the minus sign from derivation
    P_world_reconstructed = (P1 + P2) / 2

    return P_world_reconstructed


# --- Simulation Setup ---
# 1. Define a true 3D point in the world
P_world_original = np.array([2.0, 3.0, 5.0])

# 2. Define Camera Intrinsics (K)
K = np.array([
    [800, 0,   320], # [fx, 0, cx]
    [0,   800, 240], # [0, fy, cy]
    [0,   0,   1]    # [0,  0,  1]
])

# 3. Define Camera Extrinsics (Poses)
# Camera 1 is at the world origin
R1 = np.eye(3)
t1 = np.array([0.0, 0.0, 0.0])
pose1 = np.hstack((R1, t1.reshape(-1, 1)))

# Camera 2 is translated along the x-axis by 10 units
R2 = np.eye(3)
t2 = np.array([10.0, 0.0, 0.0])
pose2 = np.hstack((R2, t2.reshape(-1, 1)))

# --- Image Formation (Projecting 3D point to 2D) ---
# Project to camera 1
P_cam1 = R1.T @ (P_world_original - t1)
p_proj1 = K @ P_cam1
p1_image = p_proj1[:2] / p_proj1[2]

# Project to camera 2
P_cam2 = R2.T @ (P_world_original - t2)
p_proj2 = K @ P_cam2
p2_image = p_proj2[:2] / p_proj2[2]

print("--- Simulation Setup ---")
print(f"Original 3D Point (World Frame): {P_world_original}")
print(f"Observed 2D point in Image 1: {np.round(p1_image, 2)}")
print(f"Observed 2D point in Image 2: {np.round(p2_image, 2)}")
print("-" * 30, "\n")


# --- Main Task: Triangulation ---
P_reconstructed = triangulate_point(p1_image, p2_image, K, pose1, pose2)

print("\n--- Final Result ---")
if P_reconstructed is not None:
    print(f"Original 3D Point:      {np.round(P_world_original, 5)}")
    print(f"Reconstructed 3D Point: {np.round(P_reconstructed, 5)}")
    print("\nConclusion: The triangulation yields the point's coordinates in the common WORLD frame.")
    print("A transformation from each camera's local frame to the world frame was a necessary prerequisite.")
    print("To get the result in a camera's frame (e.g., Cam1), a final transformation would be needed: P_cam1 = R1.T @ (P_reconstructed - t1).")
    print("Therefore, the method does not *directly* yield a solution in the camera's reference frame.")