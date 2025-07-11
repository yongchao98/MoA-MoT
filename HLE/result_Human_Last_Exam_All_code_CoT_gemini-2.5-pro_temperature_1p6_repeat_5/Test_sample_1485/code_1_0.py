import numpy as np

def solve():
    """
    This function generates the rotated view of a torus based on the problem description.
    """
    # Rotation angles in degrees, as specified in the problem.
    angle_x_deg = 140
    angle_y_deg = 75
    angle_z_deg = 35
    
    # Announce the parameters being used.
    print(f"Applying rotations (X, Y, Z) in degrees: ({angle_x_deg}, {angle_y_deg}, {angle_z_deg})")
    print("-" * 50)

    def rotation_matrix_cw(axis, angle_deg):
        """Helper function to create a clockwise rotation matrix."""
        angle_rad = np.deg2rad(angle_deg)
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        if axis == 'x':
            return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
        elif axis == 'y':
            return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])
        elif axis == 'z':
            return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

    # Torus and rendering parameters
    R, r = 3.0, 1.2  # Major and minor radii
    width, height = 50, 22
    scale = 5.0
    
    # Combined rotation matrix M = Rz @ Ry @ Rx
    M = rotation_matrix_cw('z', angle_z_deg) @ rotation_matrix_cw('y', angle_y_deg) @ rotation_matrix_cw('x', angle_x_deg)

    # Screen and Z-buffer initialization
    screen = np.full((height, width), ' ')
    z_buffer = np.full((height, width), np.inf)
    chars = ['░', '▒', '▓', '█'] # Shading characters from nearest to farthest

    # First pass: Generate points and find the full range of z-values for normalization
    points_to_render = []
    z_values_for_shading = []
    u_step = 0.04
    v_step = 0.04
    
    for u in np.arange(0, 2 * np.pi, u_step):
        for v in np.arange(0, 2 * np.pi, v_step):
            # Using the special parameterization that produces the symmetric initial image
            x = (R + r * np.cos(u)) * np.cos(v)
            y = (R + r * np.cos(u)) * np.sin(v)
            z = -r * np.cos(u)  # z depends on radius, giving symmetric radial shading
            
            P_initial = np.array([x, y, z])
            P_rotated = M @ P_initial
            points_to_render.append(P_rotated)
            z_values_for_shading.append(P_rotated[2])

    # Determine z-range for shading
    if not z_values_for_shading:
      z_min, z_max = -1, 1
    else:
      z_min, z_max = min(z_values_for_shading), max(z_values_for_shading)
      
    z_range = z_max - z_min
    if z_range == 0:
        z_range = 1 # Avoid division by zero

    # Second pass: Render the points onto the screen
    for p in points_to_render:
        x_rot, y_rot, z_rot = p[0], p[1], p[2]
        
        # Project to 2D screen coordinates
        screen_x = int(width / 2 + scale * x_rot)
        screen_y = int(height / 2 - scale * y_rot)  # Y-axis is bottom-to-top

        # Z-buffering and drawing
        if 0 <= screen_x < width and 0 <= screen_y < height:
            if z_rot < z_buffer[screen_y, screen_x]:
                z_buffer[screen_y, screen_x] = z_rot
                
                # Normalize z-value and select character
                normalized_z = (z_rot - z_min) / z_range
                char_idx = min(len(chars) - 1, int(normalized_z * len(chars)))
                screen[screen_y, screen_x] = chars[char_idx]

    # Print the final result
    for row in screen:
        print("".join(row))

solve()
<<<B>>>