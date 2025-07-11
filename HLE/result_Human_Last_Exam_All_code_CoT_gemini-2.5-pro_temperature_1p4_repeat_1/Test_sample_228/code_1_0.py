import numpy as np

def triangulate_point(line1_point, line1_dir, line2_point, line2_dir):
    """
    Finds the midpoint of the shortest segment between two 3D lines.
    This point is the triangulated 3D point.
    All inputs must be in the same (world) coordinate system.
    """
    # Normalize direction vectors
    line1_dir = line1_dir / np.linalg.norm(line1_dir)
    line2_dir = line2_dir / np.linalg.norm(line2_dir)

    # Calculate the shortest line segment between the two lines
    # Based on the formula for the closest point of approach of two skew lines
    n = np.cross(line1_dir, line2_dir)
    n_mag_sq = np.dot(n, n)
    
    # Handle parallel lines case
    if np.isclose(n_mag_sq, 0):
        # For this demonstration, if lines are parallel, we can't triangulate.
        # In a real system, you'd handle this case, e.g., by returning an error.
        print("Lines are parallel, triangulation is ambiguous.")
        # As a simplification, just return the average of the starting points.
        return (line1_point + line2_point) / 2.0

    p_diff = line2_point - line1_point
    
    t1 = np.dot(np.cross(p_diff, line2_dir), n) / n_mag_sq
    t2 = np.dot(np.cross(p_diff, line1_dir), n) / n_mag_sq

    # Points on each line that are closest to the other line
    closest_point_on_line1 = line1_point + t1 * line1_dir
    closest_point_on_line2 = line2_point + t2 * line2_dir
    
    # The triangulated point is the midpoint of this shortest segment
    triangulated_point = (closest_point_on_line1 + closest_point_on_line2) / 2.0
    return triangulated_point

# --- Setup ---
# 1. Define a true 3D point in the WORLD frame.
P_world = np.array([5.0, 3.0, 20.0])

# 2. Define two camera poses in the WORLD frame.
# Camera 1 is at the world origin, looking along Z.
C1_world = np.array([0.0, 0.0, 0.0])
R1_world = np.identity(3) # No rotation relative to world

# Camera 2 is translated along the X-axis, looking along Z.
C2_world = np.array([10.0, 0.0, 0.0])
R2_world = np.identity(3) # No rotation relative to world

print(f"Original 3D Point in World Frame: {P_world}")
print("-" * 50)

# --- Simulation ---
# 3. Define the two rays (lines) in the WORLD frame.
# A ray is defined by the camera center and a direction vector pointing to the 3D point.
# These lines are what would be represented by Plucker coordinates.
ray1_origin = C1_world
ray1_direction = P_world - C1_world

ray2_origin = C2_world
ray2_direction = P_world - C2_world

# 4. Perform triangulation. The calculation uses entities defined in the WORLD frame.
P_triangulated_world = triangulate_point(ray1_origin, ray1_direction, ray2_origin, ray2_direction)

print(f"Triangulation Result: {P_triangulated_world}")
print("Note: The triangulation result is in the WORLD frame, matching the original P_world.")
print("-" * 50)

# --- Transformation ---
# 5. To get the point in a CAMERA's reference frame, a transformation is needed.
# Let's find the coordinates of P_world in Camera 2's frame.
# The transformation from world to camera is: P_cam = R_inv * (P_world - C_world)
# Since R is an orthonormal rotation matrix, R_inv = R_transpose.
R2_inv = R2_world.T
P_in_cam2_frame = R2_inv @ (P_world - C2_world)

print(f"Original 3D point's coordinates in Camera 2's frame:")
print(f"P_cam2 = R2.T @ (P_world - C2_world) = {P_in_cam2_frame}")
print("\nConclusion: The triangulation yields a point in the world frame.")
print("A separate transformation is required to express this point in a camera's reference frame.")
print("The direct result from triangulation is not in the camera frame (unless the camera is the world origin).")
