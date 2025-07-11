import numpy as np

def get_plucker_coords(p1, p2):
    """
    Computes the Plucker coordinates for a line passing through points p1 and p2.
    """
    direction = p2 - p1
    direction /= np.linalg.norm(direction) # Normalize the direction vector
    moment = np.cross(p1, p2)
    return direction, moment

def point_to_camera_frame(point_world, camera_pos_world, camera_rot_world):
    """
    Transforms a point from world coordinates to camera coordinates.
    """
    # The transformation is P_cam = R_transpose * (P_world - C_world)
    point_camera = camera_rot_world.T @ (point_world - camera_pos_world)
    return point_camera

# --- 1. Setup Scene in World Reference Frame ---
# A 3D point in the world
P_world = np.array([2.0, 3.0, 5.0])

# Camera 1: At the origin, looking along Z
C1_world = np.array([0.0, 0.0, 0.0])
R1_world = np.identity(3) # Identity matrix for rotation

# Camera 2: Shifted and rotated
C2_world = np.array([4.0, 0.0, 0.0]) # Shifted along X-axis
# Rotated -30 degrees around Y-axis
theta = np.deg2rad(-30)
R2_world = np.array([
    [np.cos(theta), 0, np.sin(theta)],
    [0, 1, 0],
    [-np.sin(theta), 0, np.cos(theta)]
])

print("--- Step-by-step Triangulation ---")
print(f"1. A 3D point exists in the world at: {P_world}")
print(f"2. Camera 1 is at {C1_world} and Camera 2 is at {C2_world} in the world frame.")

# --- 2. Define Rays in World Frame and get Plucker Coords ---
# The ray for each camera is the line from its center to the 3D point.
# These lines *must* be defined in the common World Frame.
print("\n3. Lines from each camera's center to the point are defined in the World Frame.")
line1_dir_world, line1_moment_world = get_plucker_coords(C1_world, P_world)
line2_dir_world, line2_moment_world = get_plucker_coords(C2_world, P_world)
print(f"   - Line 1 (Plucker direction) in world frame: {np.round(line1_dir_world, 2)}")
print(f"   - Line 2 (Plucker direction) in world frame: {np.round(line2_dir_world, 2)}")

# --- 3. Triangulation (Conceptual) ---
# A real triangulation algorithm finds the point of closest approach between two skew lines.
# For this ideal example, the intersection is simply our original point, P_world.
# The crucial concept is that the triangulation algorithm operates on the world-frame
# line representations, so its output is also in the world frame.
triangulated_point_world = P_world # In a real case, this would be computed
print(f"\n4. Triangulation finds the intersection of the lines.")
print(f"   The DIRECT result of triangulation is the point in the World Frame.")
print(f"   Triangulated Point (World Coords): [{triangulated_point_world[0]}, {triangulated_point_world[1]}, {triangulated_point_world[2]}]")

# --- 4. Transform to Camera Frame ---
# To get the point's coordinates in Camera 2's reference frame, we need to apply
# the transformation defined by Camera 2's pose (extrinsic parameters).
P_cam2 = point_to_camera_frame(triangulated_point_world, C2_world, R2_world)

print(f"\n5. To get the point's coordinates in a specific camera's frame, a transformation is needed.")
print(f"   Applying the transformation for Camera 2 (Rotation and Translation)...")
print(f"   Point in Camera 2 Coords: [{P_cam2[0]:.2f}, {P_cam2[1]:.2f}, {P_cam2[2]:.2f}]")

print("\nConclusion: Triangulation directly yields a solution in the common World Frame.")
print("A separate transformation is required to represent that solution in a camera's local frame.")
