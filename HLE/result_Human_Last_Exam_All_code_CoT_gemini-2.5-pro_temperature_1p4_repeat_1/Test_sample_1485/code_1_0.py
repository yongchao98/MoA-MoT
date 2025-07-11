import numpy as np

def generate_rotated_torus_view():
    """
    This function generates and prints an ASCII representation of a torus
    after being rotated by specific angles around the X, Y, and Z axes.
    """
    # --- Parameters ---
    # Screen and torus parameters
    width, height = 50, 22
    major_radius, minor_radius = 8, 3.5
    scale_factor = 3.5

    # Rotation angles in degrees, as specified in the problem
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35

    # --- Setup ---
    # Convert angles to radians for trigonometric functions
    theta_x = np.radians(rot_x_deg)
    theta_y = np.radians(rot_y_deg)
    theta_z = np.radians(rot_z_deg)

    # Define clockwise rotation matrices
    rot_matrix_x = np.array([
        [1, 0, 0],
        [0, np.cos(theta_x), np.sin(theta_x)],
        [0, -np.sin(theta_x), np.cos(theta_x)]
    ])
    rot_matrix_y = np.array([
        [np.cos(theta_y), 0, -np.sin(theta_y)],
        [0, 1, 0],
        [np.sin(theta_y), 0, np.cos(theta_y)]
    ])
    rot_matrix_z = np.array([
        [np.cos(theta_z), np.sin(theta_z), 0],
        [-np.sin(theta_z), np.cos(theta_z), 0],
        [0, 0, 1]
    ])
    
    # Combine rotations: applied in order X, then Y, then Z
    combined_rotation_matrix = rot_matrix_z @ rot_matrix_y @ rot_matrix_x

    # Initialize buffers for rendering
    z_buffer = np.full((height, width), np.inf)
    depth_map = np.full((height, width), np.nan)

    # --- 3D Calculation and Projection ---
    # Generate points on the torus surface, rotate them, and project to screen
    u_step, v_step = 0.05, 0.02
    for u in np.arange(0, 2 * np.pi, u_step):
        for v in np.arange(0, 2 * np.pi, v_step):
            # Initial point for a torus in the YZ plane (hole along X-axis)
            x0 = minor_radius * np.sin(u)
            y0 = (major_radius + minor_radius * np.cos(u)) * np.cos(v)
            z0 = (major_radius + minor_radius * np.cos(u)) * np.sin(v)
            initial_point = np.array([x0, y0, z0])
            
            # Rotate the point
            rotated_point = combined_rotation_matrix @ initial_point
            
            # Project onto screen coordinates
            x_proj, y_proj, z_depth = rotated_point[0], rotated_point[1], rotated_point[2]
            screen_x = int(width / 2 + scale_factor * x_proj)
            screen_y = int(height / 2 - scale_factor * y_proj) # Y is inverted on screens

            # Z-buffering to handle occlusion
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if z_depth < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = z_depth
                    depth_map[screen_y, screen_x] = z_depth

    # --- Shading and Rendering ---
    # Determine depth range for shading
    visible_depths = depth_map[~np.isnan(depth_map)]
    if visible_depths.size == 0:
        min_z, z_range = 0, 1
    else:
        min_z = np.min(visible_depths)
        max_z = np.max(visible_depths)
        z_range = max_z - min_z
        if z_range == 0:
            z_range = 1 # Avoid division by zero

    # Define shades: from nearest (░) to farthest (█)
    shades = ['░', '▒', '▓', '█']
    output_buffer = np.full((height, width), ' ')

    # Populate the output buffer with shaded characters
    for r_idx in range(height):
        for c_idx in range(width):
            depth = depth_map[r_idx, c_idx]
            if not np.isnan(depth):
                normalized_depth = (depth - min_z) / z_range
                shade_index = int(normalized_depth * 3.99)
                output_buffer[r_idx, c_idx] = shades[shade_index]

    # --- Print the Result ---
    print(f"The view after rotation by (X={rot_x_deg}, Y={rot_y_deg}, Z={rot_z_deg}) degrees is:")
    print("".join([' ' * (width + 4)])) # Top padding
    for row in output_buffer:
        print("  " + "".join(row) + "  ")
    print("".join([' ' * (width + 4)])) # Bottom padding

if __name__ == '__main__':
    generate_rotated_torus_view()