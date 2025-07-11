import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    """
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def triangulate_point(p1, p2, P1, P2):
    """
    Triangulates a 3D point from two 2D points and two projection matrices.
    This is a classic direct linear transform (DLT) approach, which serves as a comparison.
    A more direct method for two lines is used in the main script.
    """
    A = np.array([
        p1[0] * P1[2, :] - P1[0, :],
        p1[1] * P1[2, :] - P1[1, :],
        p2[0] * P2[2, :] - P2[0, :],
        p2[1] * P2[2, :] - P2[1, :]
    ])
    A = A.reshape((4, 4))
    _, _, V = np.linalg.svd(A)
    X = V[-1, :4]
    return X / X[3]


def main():
    """
    Main function to simulate triangulation and demonstrate its limitations.
    """
    # 1. Simulate a 3D scene and cameras
    # Intrinsics matrix (assuming same for both cameras)
    K = np.array([[800, 0, 320], [0, 800, 240], [0, 0, 1]])
    
    # Extrinsics for Camera 1 (at the world origin)
    R1 = np.identity(3)
    t1 = np.zeros((3, 1))
    
    # Extrinsics for Camera 2
    # Rotate around Y axis and translate along X
    angle = np.deg2rad(15)
    R2 = np.array([
        [np.cos(angle), 0, np.sin(angle)],
        [0, 1, 0],
        [-np.sin(angle), 0, np.cos(angle)]
    ])
    t2 = np.array([[-2], [0], [0.5]])

    # Define the true 3D point in the world
    P_world_true = np.array([[1.5], [0.5], [4.0]])
    print(f"Original 3D Point:\n{P_world_true.flatten()}\n")

    # 2. Project the point onto camera image planes
    # Projection matrices
    P1 = K @ np.hstack((R1, t1))
    P2 = K @ np.hstack((R2, t2))

    # Project P_world into each camera's view
    p1_homogeneous = P1 @ np.vstack((P_world_true, [1]))
    p1_img = (p1_homogeneous / p1_homogeneous[2])[:2].flatten()
    
    p2_homogeneous = P2 @ np.vstack((P_world_true, [1]))
    p2_img = (p2_homogeneous / p2_homogeneous[2])[:2].flatten()
    
    print(f"Perfect 2D projection in Cam1: {p1_img}")
    print(f"Perfect 2D projection in Cam2: {p2_img}\n")

    # 3. Introduce noise to simulate measurement error
    noise_level = 0.5  # pixels
    p1_noisy = p1_img + np.random.uniform(-noise_level, noise_level, 2)
    p2_noisy = p2_img + np.random.uniform(-noise_level, noise_level, 2)
    
    print(f"Noisy 2D projection in Cam1: {p1_noisy}")
    print(f"Noisy 2D projection in Cam2: {p2_noisy}\n")

    # 4. Reconstruct 3D rays from the noisy 2D points
    # A ray is defined by its origin (camera center) and direction.
    
    # Camera centers in world coordinates
    C1_world = -R1.T @ t1
    C2_world = -R2.T @ t2
    
    # Back-project noisy 2D points to get ray directions in world coordinates
    d1_world = R1.T @ np.linalg.inv(K) @ np.array([p1_noisy[0], p1_noisy[1], 1])
    d1_world = d1_world / np.linalg.norm(d1_world) # Normalize
    
    d2_world = R2.T @ np.linalg.inv(K) @ np.array([p2_noisy[0], p2_noisy[1], 1])
    d2_world = d2_world / np.linalg.norm(d2_world) # Normalize

    # These rays are L1(s) = C1 + s*d1 and L2(t) = C2 + t*d2

    # Representing lines using Plücker coordinates (implicitly) to check for intersection
    # Plücker coordinates are L = (direction, moment). Moment m = origin x direction.
    m1_world = np.cross(C1_world.flatten(), d1_world)
    m2_world = np.cross(C2_world.flatten(), d2_world)

    print("Reconstructed 3D rays (in world coordinates):")
    print(f"Ray 1 Origin (C1): {C1_world.flatten()}")
    print(f"Ray 1 Direction (d1): {d1_world.flatten()}")
    print(f"Ray 2 Origin (C2): {C2_world.flatten()}")
    print(f"Ray 2 Direction (d2): {d2_world.flatten()}\n")

    # 5. Check for intersection using the side operator from Plücker geometry
    # Two lines L1=(d1,m1) and L2=(d2,m2) intersect if and only if d1·m2 + d2·m1 = 0.
    # This value represents the signed shortest distance between the lines.
    intersection_test = d1_world @ m2_world + d2_world @ m1_world
    
    print("--- Intersection Test ---")
    print(f"Plücker line intersection condition (d1·m2 + d2·m1): {intersection_test:.6f}")
    if np.isclose(intersection_test, 0):
        print("Result: The lines appear to intersect perfectly (this is rare with noise).")
    else:
        print("Result: The lines are skew and DO NOT intersect. A direct solution is not possible.\n")

    # 6. Triangulate the point that is closest to both non-intersecting lines
    # Find the midpoint of the shortest line segment connecting the two rays.
    # We solve the system [d1, -d2][s; t] = C2 - C1 for the closest approach,
    # which can be done via least squares.
    
    A = np.vstack((d1_world, -d2_world)).T
    b = (C2_world - C1_world).flatten()
    
    # Least squares solution for parameters s and t
    # (A.T @ A) @ x = A.T @ b
    try:
        st_params, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        s, t = st_params
        
        # Points on each line closest to the other line
        P1_closest = C1_world.flatten() + s * d1_world
        P2_closest = C2_world.flatten() + t * d2_world
        
        # The triangulated point is the midpoint of this shortest segment
        P_world_triangulated = (P1_closest + P2_closest) / 2

        print("--- Triangulation of Non-Intersecting Lines ---")
        print(f"Closest point on Ray 1: {P1_closest}")
        print(f"Closest point on Ray 2: {P2_closest}")
        print(f"Final Triangulated Point (Midpoint): {P_world_triangulated}")

        error = np.linalg.norm(P_world_triangulated - P_world_true.flatten())
        print(f"\nReconstruction Error: {error:.4f}")

    except np.linalg.LinAlgError:
        print("Triangulation failed. This can happen if the rays are parallel.")


if __name__ == '__main__':
    main()