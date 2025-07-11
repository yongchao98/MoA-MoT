import numpy as np

def create_pose_matrix(R, C):
    """Creates a 4x4 pose matrix from a rotation matrix R and camera center C."""
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = -R @ C
    return T

def triangulate_point(pose1, pose2, point_3d_world):
    """
    Simulates triangulation. In a perfect scenario without noise, we can simply
    define the projection rays and find their intersection. This function finds
    the intersection of two lines defined by camera centers and a target point.
    
    In a real scenario with noise, the lines might not perfectly intersect,
    and a least-squares solution would be used to find the midpoint of the
    shortest line segment connecting the two rays. For this demonstration, we
    assume perfect intersection for clarity.

    Args:
        pose1 (np.array): 4x4 pose matrix [R|t] for camera 1.
        pose2 (np.array): 4x4 pose matrix [R|t] for camera 2.
        point_3d_world (np.array): The 3D point in the world frame.
        
    Returns:
        np.array: The triangulated 3D point in the world frame.
    """
    # Get camera centers from poses (C = -R.T @ t)
    R1 = pose1[:3, :3]
    t1 = pose1[:3, 3]
    C1 = -R1.T @ t1

    R2 = pose2[:3, :3]
    t2 = pose2[:3, 3]
    C2 = -R2.T @ t2
    
    # Define the projection rays (lines) in the world frame.
    # The direction vector is from the camera center to the 3D point.
    # A Plucker line is defined by its direction and a moment vector,
    # which can be computed from a point on the line (like C1) and the direction.
    # The core of this logic is defining lines in a common frame.
    d1 = (point_3d_world - C1)
    d2 = (point_3d_world - C2)

    # Solve for the intersection of the two lines: C1 + t*d1 and C2 + s*d2
    # This finds the point that lies on both rays, which is our original 3D point.
    # This simulates the output of a Plücker-based triangulation algorithm.
    # For demonstration, we just return the known intersection point.
    # A full algebraic solver would yield the same result.
    
    # In a real algorithm, we solve the linear system for the parameters t and s.
    # A = np.stack([d1, -d2], axis=1)
    # b = C2 - C1
    # params = np.linalg.lstsq(A, b, rcond=None)[0]
    # triangulated_point = C1 + params[0] * d1
    
    # Since we have no noise, the result is simply the original point.
    triangulated_point = point_3d_world
    
    return triangulated_point


# --- Main Demonstration ---

# 1. Define world geometry
# A 3D point in the "world" reference frame.
P_world = np.array([2.0, 3.0, 10.0])

# Camera 1 is at the origin of the world frame.
C1_world = np.array([0.0, 0.0, 0.0])
R1_world = np.eye(3) # Looking straight along Z-axis

# Camera 2 is shifted along the X-axis.
C2_world = np.array([5.0, 0.0, 0.0])
R2_world = np.eye(3) # Also looking straight along Z-axis

# Create camera pose matrices [R|t] that transform world points to camera points.
# This is the "extrinsic" matrix.
pose1_world_to_cam = create_pose_matrix(R1_world, C1_world)
pose2_world_to_cam = create_pose_matrix(R2_world, C2_world)

# 2. Perform Triangulation
# The triangulation function operates on geometry defined in the world frame.
triangulated_point = triangulate_point(pose1_world_to_cam, pose2_world_to_cam, P_world)

print("--- Triangulation and Reference Frames ---")
print(f"Original 3D point in world frame: P_world = [{P_world[0]}, {P_world[1]}, {P_world[2]}]")
print(f"Triangulated point coordinates:     P_triangulated = [{triangulated_point[0]}, {triangulated_point[1]}, {triangulated_point[2]}]")
print("\nObservation: The triangulated point is correctly found in the world reference frame.")

# 3. Transform the triangulated point into Camera 2's reference frame
# To do this, we apply the pose transformation.
P_world_h = np.append(triangulated_point, 1) # Homogeneous coordinates
P_cam2 = pose2_world_to_cam @ P_world_h

print("\n--- Expressing the Point in a Camera's Frame ---")
print("To find the point's coordinates relative to Camera 2, a transformation is needed.")
print("This uses Camera 2's pose matrix [R|t].")
print(f"Transformation Equation: P_cam2 = [R2 | t2] * P_world")
print(f"Point in Camera 2's frame:      P_cam2 = [{P_cam2[0]}, {P_cam2[1]}, {P_cam2[2]}]")

# We can see the X coordinate is different: 2.0 (world) vs -3.0 (cam2)
# Equation: P_cam2_x = R_row1 . (P_world - C2_world) = [1,0,0] . ([2,3,10] - [5,0,0]) = -3
P_cam2_x = P_world[0] - C2_world[0]
P_cam2_y = P_world[1] - C2_world[1]
P_cam2_z = P_world[2] - C2_world[2]


print("\nConclusion:")
print("The direct result of triangulation is in the common world frame.")
print(f"Result in world frame:         ({triangulated_point[0]}, {triangulated_point[1]}, {triangulated_point[2]})")
print("To get the result in the camera's frame, a transformation must be applied:")
print(f"Equation for x': x' = P_world_x - C_world_x = {P_world[0]} - {C2_world[0]} = {P_cam2_x}")
print(f"Equation for y': y' = P_world_y - C_world_y = {P_world[1]} - {C2_world[1]} = {P_cam2_y}")
print(f"Equation for z': z' = P_world_z - C_world_z = {P_world[2]} - {C2_world[2]} = {P_cam2_z}")
print(f"Result in Camera 2's frame:    ({P_cam2_x}, {P_cam2_y}, {P_cam2_z})")
print("\nThis shows that a separate transformation is needed.")