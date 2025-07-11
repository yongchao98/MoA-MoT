import numpy as np

def triangulate_point(pt_3d_original, P1, P2):
    """
    Simulates triangulation.
    1. Projects a 3D point into two cameras.
    2. Uses the 2D projections and camera matrices to triangulate the 3D point.
    """
    # Ensure the original point is in homogeneous coordinates for projection
    pt_3d_homogeneous = np.append(pt_3d_original, 1)

    # Project the 3D point onto both camera image planes
    pt_2d_cam1 = P1 @ pt_3d_homogeneous
    pt_2d_cam2 = P2 @ pt_3d_homogeneous

    # Normalize the 2D points (divide by the last coordinate)
    pt_2d_cam1 = pt_2d_cam1 / pt_2d_cam1[2]
    pt_2d_cam2 = pt_2d_cam2 / pt_2d_cam2[2]

    print("--- Simulation Inputs ---")
    print(f"Original 3D point (in Camera 1's frame): {pt_3d_original}")
    print(f"Projection in Camera 1: {pt_2d_cam1[:2].flatten()}")
    print(f"Projection in Camera 2: {pt_2d_cam2[:2].flatten()}\n")

    # Now, use the camera matrices and the 2D points to triangulate
    # This simulates having only the 2D measurements and camera poses
    # Note: cv2.triangulatePoints expects 2xN or 2x1 points
    pt_2d_cam1_for_triangulation = pt_2d_cam1[:2].reshape(2, 1)
    pt_2d_cam2_for_triangulation = pt_2d_cam2[:2].reshape(2, 1)

    # In a real system using Plucker coordinates, this step would be replaced
    # by forming the linear system from the Plucker line incidence constraints
    # and solving for the 3D point X. Here we use a standard library function
    # which solves the same underlying geometric problem.
    # The output is in the same coordinate frame that P1 and P2 are defined in.
    pt_3d_triangulated_hom = cv2_triangulate_points(P1, P2, pt_2d_cam1_for_triangulation, pt_2d_cam2_for_triangulation)

    # Convert from homogeneous to Cartesian coordinates
    pt_3d_triangulated = pt_3d_triangulated_hom[:3] / pt_3d_triangulated_hom[3]
    
    return pt_3d_triangulated.flatten()

def cv2_triangulate_points(P1, P2, points1, points2):
    """
    A simplified wrapper for a triangulation function. In a real scenario,
    one might use an OpenCV implementation. Here we recreate the linear algebra
    for clarity, which is analogous to the Plucker coordinate method.
    For a point x = PX, we can write this as x_i = P_i^T X where P_i^T is the
    i-th row of P. Then cross product (x_i X P_i^T X) = 0.
    This gives two linear equations in X per view. We stack them to solve AX=0.
    """
    A = np.zeros((4, 4))
    A[0] = points1[0] * P1[2] - P1[0]
    A[1] = points1[1] * P1[2] - P1[1]
    A[2] = points2[0] * P2[2] - P2[0]
    A[3] = points2[1] * P2[2] - P2[1]
    
    # Solve AX=0 using Singular Value Decomposition
    _u, _s, vh = np.linalg.svd(A)
    # The solution is the last column of V (or last row of V^T)
    X = vh[-1]
    return X


# --- Main execution ---
# 1. Define Camera Intrinsics (e.g., a simple camera)
fx, fy = 600, 600
cx, cy = 320, 240
K = np.array([[fx, 0, cx],
              [0, fy, cy],
              [0, 0, 1]])

# 2. Define Camera Poses.
# We set the "world" coordinate system to be the same as Camera 1's system.
# Pose of Camera 1 (relative to the world/itself) is identity.
R1 = np.eye(3)
t1 = np.zeros((3, 1))

# Pose of Camera 2 (relative to Camera 1) - translated 0.5 units along X-axis
R2 = np.eye(3)
t2 = np.array([[0.5], [0], [0]]) # 0.5m baseline

# 3. Create Projection Matrices
# Projection Matrix P = K * [R | t]
P1 = K @ np.hstack((R1, t1)) # This is K @ [I | 0]
P2 = K @ np.hstack((R2, t2))

# 4. Define a sample 3D point IN THE COORDINATE FRAME OF CAMERA 1.
# Since we aligned the world with Camera 1, this is also its world coordinate.
pt_3d_original = np.array([1.0, 0.5, 5.0]) # x=1.0, y=0.5, z=5.0 meters away

# 5. Run the simulation
pt_3d_reconstructed = triangulate_point(pt_3d_original, P1, P2)

# 6. Print and Compare Results
print("--- Triangulation Result ---")
print("This result is directly in the coordinate frame of Camera 1 because")
print("we defined the triangulation problem with Camera 1 as the reference origin.\n")

print(f"Original 3D point:      x={pt_3d_original[0]:.6f}, y={pt_3d_original[1]:.6f}, z={pt_3d_original[2]:.6f}")
print(f"Reconstructed 3D point: x={pt_3d_reconstructed[0]:.6f}, y={pt_3d_reconstructed[1]:.6f}, z={pt_3d_reconstructed[2]:.6f}")

<<<C>>>