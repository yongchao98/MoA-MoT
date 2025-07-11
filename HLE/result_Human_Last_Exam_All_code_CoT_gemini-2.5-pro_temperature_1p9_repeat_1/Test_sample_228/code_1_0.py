import numpy as np

def skew(v):
    """Converts a 3-element vector to its skew-symmetric matrix form."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def plucker_from_points(A, B):
    """Creates Plucker coordinates for a line passing through points A and B."""
    A = np.asarray(A)
    B = np.asarray(B)
    direction = B - A
    moment = np.cross(A, B)
    # Normalize for consistency
    norm_dir = np.linalg.norm(direction)
    if norm_dir == 0:
        return np.zeros(6)
    return np.concatenate((direction / norm_dir, moment / norm_dir))

def main():
    """
    Demonstrates the necessity of a coordinate transformation for triangulation
    using Plucker coordinates.
    """
    print("Step 1: Define scene geometry (world, cameras, 3D point).")
    # Camera 1 is at the origin, its frame is the World Frame.
    C1_world = np.array([0., 0., 0.])
    
    # Camera 2 is translated and rotated relative to the world.
    # Pose of Camera 2 in the world (transforms points from cam2 frame to world frame)
    theta = np.pi / 4  # 45-degree rotation around Y axis
    R_cam2_to_world = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])
    t_cam2_in_world = np.array([10., 1., 2.]) # Position of Cam2 center in world

    # Define a single 3D point in the World Frame.
    P_world = np.array([2., 3., 7.])
    print(f"Original 3D point in World Frame: {np.round(P_world, 2)}\n")

    print("Step 2: Define Ray 1 in the World Frame.")
    # Ray 1 is from the center of Cam1 (origin) through P_world.
    L1_world = plucker_from_points(C1_world, P_world)
    print(f"Ray 1 (from Cam1) in World Frame coordinates: L1_world = [{', '.join([f'{x:.2f}' for x in L1_world])}]")

    print("\nStep 3: Define Ray 2 in its local camera frame.")
    # To get Ray 2 in its own frame, we need P's coordinates relative to Cam2.
    # First, transform P from world coordinates to Cam2's local coordinates.
    P_in_cam2_frame = R_cam2_to_world.T @ (P_world - t_cam2_in_world)
    
    # In its own frame, the center of Cam2 is the origin.
    C2_local = np.array([0., 0., 0.])
    L2_local = plucker_from_points(C2_local, P_in_cam2_frame)
    print(f"Ray 2 (from Cam2) in its LOCAL camera frame: L2_local = [{', '.join([f'{x:.2f}' for x in L2_local])}]")
    print("Notice L1_world and L2_local have different values as they are in different coordinate systems.\n")
    
    print("Step 4: Show the limitation. We cannot combine L1_world and L2_local.")
    print("A transformation is required to bring L2 into the same frame as L1.\n")
    
    print("Step 5: Apply the Plucker line transformation.")
    # The transformation matrix for Plucker lines from cam2 frame to world frame.
    plucker_transform = np.zeros((6, 6))
    plucker_transform[0:3, 0:3] = R_cam2_to_world
    plucker_transform[3:6, 3:6] = R_cam2_to_world
    plucker_transform[3:6, 0:3] = skew(t_cam2_in_world) @ R_cam2_to_world

    # Apply the transformation to the local ray coordinates.
    L2_world = plucker_transform @ L2_local
    print(f"Ray 2 coordinates after transformation to World Frame:")
    
    # Final step of the triangulation would use L1_world and L2_world.
    # The "equation" for triangulation is built from these two commensurate vectors.
    print(f"Final Triangulation Input Equation:")
    L1_parts = [f"{x:.2f}" for x in L1_world]
    L2_parts = [f"{x:.2f}" for x in L2_world]
    print(f"   Ray 1 (World): [{', '.join(L1_parts)}]")
    print(f" + Ray 2 (World): [{', '.join(L2_parts)}]")
    print("-------------------------------------------------------------------------")
    print("Conclusion: A transformation is a required step before triangulation can proceed.")


if __name__ == '__main__':
    main()
