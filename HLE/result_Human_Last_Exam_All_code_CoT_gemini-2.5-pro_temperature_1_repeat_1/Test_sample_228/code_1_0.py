import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix.
    """
    if len(v) != 3:
        raise ValueError("Input vector must have 3 elements")
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def main():
    """
    Demonstrates 3D point triangulation using Plucker coordinates,
    highlighting the need for a coordinate transformation.
    """
    print("--- 3D Triangulation with Plücker Coordinates ---")
    
    # 1. Define the scene setup (a 3D point and two cameras)
    # The world coordinate system is aligned with Camera 1's frame.
    
    # The 3D point in the world frame
    P_world = np.array([2.0, 3.0, 10.0])
    print(f"Original 3D Point (World Frame): {P_world}")

    # Camera 1 is at the origin of the world frame
    C1 = np.array([0.0, 0.0, 0.0])
    R1 = np.identity(3)

    # Camera 2 is translated and rotated relative to Camera 1 (the world frame)
    # This defines the transformation from Camera 2's frame to Camera 1's frame
    C2 = np.array([1.0, 0.0, 0.0]) # Translation vector
    # A 15-degree rotation around the Y-axis
    theta = np.deg2rad(15)
    R2 = np.array([[np.cos(theta), 0, np.sin(theta)],
                   [0, 1, 0],
                   [-np.sin(theta), 0, np.cos(theta)]])

    # 2. Define the rays (as Plücker lines) in their respective camera frames.
    # A ray from the camera center to the 3D point is defined by:
    # Direction (l): normalized vector from camera center to the point
    # Moment (m): cross product of a point on the line (camera center) and the direction.
    # Since the camera center is the origin of its own frame, m = 0 x l = 0.
    
    # Ray 1 in Camera 1's frame
    l1_cam1 = (P_world - C1) / np.linalg.norm(P_world - C1)
    m1_cam1 = np.cross(C1, l1_cam1) # Will be [0,0,0]
    L1_cam1 = np.hstack((l1_cam1, m1_cam1))
    print(f"\nRay 1 (in Cam1 Frame) L1 = [l1, m1]: {np.round(L1_cam1, 3)}")
    
    # To get Ray 2, first express the 3D point in Camera 2's frame
    P_cam2 = R2.T @ (P_world - C2)
    l2_cam2 = P_cam2 / np.linalg.norm(P_cam2)
    m2_cam2 = np.array([0.0, 0.0, 0.0]) # Moment is zero in its own frame
    L2_cam2 = np.hstack((l2_cam2, m2_cam2))
    print(f"Ray 2 (in Cam2 Frame) L2 = [l2, m2]: {np.round(L2_cam2, 3)}")

    # 3. THE CRITICAL STEP: Transform Ray 2 into Camera 1's frame
    # A direct intersection of L1_cam1 and L2_cam2 is not possible
    # as they exist in different coordinate systems.
    print("\nApplying transformation to bring Ray 2 into Camera 1's coordinate frame...")
    l2_in_cam1_frame = R2 @ l2_cam2
    # The moment transformation rule: m_new = R @ m_old + skew(t) @ R @ l_old
    m2_in_cam1_frame = R2 @ m2_cam2 + skew(C2) @ (R2 @ l2_cam2)
    
    L2_transformed = np.hstack((l2_in_cam1_frame, m2_in_cam1_frame))
    print(f"Transformed Ray 2 (in Cam1 Frame) L2' = [l2', m2']: {np.round(L2_transformed, 3)}")
    
    # 4. Triangulate the point using the two lines now in the SAME frame
    # The intersection point P satisfies: m1 = P x l1 and m2 = P x l2
    # This can be rewritten as: skew(l1)P = m1 and skew(l2)P = m2
    # We solve this overdetermined system using least squares.
    
    A = np.vstack((skew(l1_cam1), skew(l2_in_cam1_frame)))
    b = np.hstack((m1_cam1, m2_in_cam1_frame))
    
    # Solve Ax = b for x (the 3D point P)
    P_triangulated, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    
    print("\n--- Triangulation Result ---")
    print(f"Original point P: {P_world[0]:.4f}, {P_world[1]:.4f}, {P_world[2]:.4f}")
    print(f"Triangulated point P': {P_triangulated[0]:.4f}, {P_triangulated[1]:.4f}, {P_triangulated[2]:.4f}")
    print("\nThe result shows that a transformation was essential to solve for the point.")

if __name__ == '__main__':
    main()
