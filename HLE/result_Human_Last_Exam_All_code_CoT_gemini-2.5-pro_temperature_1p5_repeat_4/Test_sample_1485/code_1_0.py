import numpy as np

def solve_torus_rotation():
    """
    This function generates a 2D ASCII representation of a torus
    after a specified 3D rotation.
    """
    # Input rotation values from the problem
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35

    # --- 1. Settings ---
    # Screen dimensions
    width, height = 52, 22
    # Torus parameters (major and minor radii)
    R, r = 2.0, 1.0
    # Scale factor for projection
    scale = 11

    # Convert degrees to radians. The problem defines positive rotation as clockwise.
    # Standard rotation matrices are counter-clockwise, so we use negative angles.
    rot_x = np.radians(-rot_x_deg)
    rot_y = np.radians(-rot_y_deg)
    rot_z = np.radians(-rot_z_deg)

    # --- 2. Rotation Matrices (for counter-clockwise rotation) ---
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rot_x), -np.sin(rot_x)],
        [0, np.sin(rot_x), np.cos(rot_x)]
    ])
    Ry = np.array([
        [np.cos(rot_y), 0, np.sin(rot_y)],
        [0, 1, 0],
        [-np.sin(rot_y), 0, np.cos(rot_y)]
    ])
    Rz = np.array([
        [np.cos(rot_z), -np.sin(rot_z), 0],
        [np.sin(rot_z), np.cos(rot_z), 0],
        [0, 0, 1]
    ])
    # Combined rotation matrix (applied as Rz * Ry * Rx)
    rotation_matrix = Rz @ Ry @ Rx

    # --- 3. Canvas and Z-buffer Initialization ---
    # The canvas will hold the final ASCII characters.
    canvas = np.full((height, width), ' ')
    # The Z-buffer stores the depth of the closest point for each pixel.
    z_buffer = np.full((height, width), np.inf)

    # --- 4. Generate, Rotate, and Project Torus Points ---
    # Iterate over the torus surface using two angles, theta and phi.
    theta_step = 0.07  # Angle for the major circle
    phi_step = 0.07    # Angle for the minor circle (tube)

    for theta in np.arange(0, 2 * np.pi, theta_step):
        for phi in np.arange(0, 2 * np.pi, phi_step):
            # Parametric equations for a point on the torus surface
            x = (R + r * np.cos(theta)) * np.cos(phi)
            y = (R + r * np.cos(theta)) * np.sin(phi)
            z = r * np.sin(theta)
            p = np.array([x, y, z])

            # Apply the combined rotation to the point
            p_rotated = rotation_matrix @ p
            x_rot, y_rot, z_rot = p_rotated[0], p_rotated[1], p_rotated[2]

            # Project the 3D point to 2D screen coordinates (orthographic projection)
            # Y-axis is inverted because screen coordinates start from top-left.
            screen_x = int(width / 2 + scale * x_rot)
            screen_y = int(height / 2 - scale * y_rot)

            # --- 5. Z-buffering (handle occlusion) ---
            if 0 <= screen_x < width and 0 <= screen_y < height:
                # If this point is closer than what's already at this pixel...
                if z_rot < z_buffer[screen_y, screen_x]:
                    # ...update the z-buffer with the new, closer depth.
                    z_buffer[screen_y, screen_x] = z_rot

    # --- 6. Determine Shading and Render to Canvas ---
    # Find min and max depth values to normalize them for shading.
    min_z = np.min(z_buffer[z_buffer != np.inf])
    max_z = np.max(z_buffer[z_buffer != np.inf])
    z_range = max_z - min_z if max_z > min_z else 1.0
    
    # Lighter shades (░) for nearer points (smaller z). Darker shades (█) for farther points.
    shades = '░▒▓█'

    for y_idx in range(height):
        for x_idx in range(width):
            if z_buffer[y_idx, x_idx] != np.inf:
                # Normalize depth to a 0-1 range.
                normalized_z = (z_buffer[y_idx, x_idx] - min_z) / z_range
                # Map the normalized depth to a character index.
                shade_index = int(normalized_z * (len(shades) - 0.001))
                canvas[y_idx, x_idx] = shades[shade_index]

    # --- 7. Print the Result ---
    print(f"The final view of the torus after a rotation of ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) around the X, Y and Z axis is:")
    print("\n" * 3) # Add spacing for better presentation
    for row in canvas:
        print("".join(row))

solve_torus_rotation()