import numpy as np

def generate_rotated_torus_view():
    """
    This function generates and prints an ASCII view of a torus
    after applying specified rotations.
    """
    # 1. SETUP: Constants for the torus, screen, and rotations
    WIDTH, HEIGHT = 50, 22
    R1, R2 = 2, 1  # Major and minor radii of the torus
    SCALE = 7      # Scale for projection size
    ASPECT_RATIO = 0.5 # Correction for character aspect ratio

    # Rotation angles in degrees from the problem
    angle_x_deg = 140
    angle_y_deg = 75
    angle_z_deg = 35

    print(f"Applying rotations of ({angle_x_deg}, {angle_y_deg}, {angle_z_deg}) degrees around the X, Y, and Z axes.")

    # Convert angles to radians for numpy functions
    ax, ay, az = np.deg2rad([angle_x_deg, angle_y_deg, angle_z_deg])

    # 2. ROTATION MATRICES: Create clockwise rotation matrices
    # Clockwise rotation around X-axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(ax), np.sin(ax)],
        [0, -np.sin(ax), np.cos(ax)]
    ])
    # Clockwise rotation around Y-axis
    Ry = np.array([
        [np.cos(ay), 0, -np.sin(ay)],
        [0, 1, 0],
        [np.sin(ay), 0, np.cos(ay)]
    ])
    # Clockwise rotation around Z-axis
    Rz = np.array([
        [np.cos(az), np.sin(az), 0],
        [-np.sin(az), np.cos(az), 0],
        [0, 0, 1]
    ])

    # Combine matrices for a single transformation. Order is Rz @ Ry @ Rx.
    R = Rz @ Ry @ Rx

    # 3. BUFFERS: Initialize screen and z-buffer
    z_buffer = np.full((HEIGHT, WIDTH), np.inf)

    # 4. TRANSFORMATION LOOP: Generate points, rotate, and project
    theta_step = 0.05
    phi_step = 0.05
    for phi in np.arange(0, 2 * np.pi, phi_step):
        for theta in np.arange(0, 2 * np.pi, theta_step):
            # Parametric equations for a torus in the XY plane
            x = (R1 + R2 * np.cos(theta)) * np.cos(phi)
            y = (R1 + R2 * np.cos(theta)) * np.sin(phi)
            z = R2 * np.sin(theta)
            
            # Apply the combined rotation to the initial point
            x_rot, y_rot, z_rot = R @ np.array([x, y, z])
            
            # Project the 3D point to 2D screen coordinates
            proj_x = int(WIDTH / 2 + SCALE * x_rot)
            proj_y = int(HEIGHT / 2 - SCALE * y_rot * ASPECT_RATIO) 
            
            # Z-buffering check
            if 0 <= proj_x < WIDTH and 0 <= proj_y < HEIGHT:
                if z_rot < z_buffer[proj_y, proj_x]:
                    z_buffer[proj_y, proj_x] = z_rot

    # 5. RENDER: Convert the z-buffer to a shaded ASCII image
    screen = np.full((HEIGHT, WIDTH), ' ')
    valid_depths = z_buffer[z_buffer != np.inf]
    
    if len(valid_depths) > 0:
        min_z = np.min(valid_depths)
        max_z = np.max(valid_depths)

        # Shades from nearest (lightest) to farthest (darkest)
        shades = ['░', '▒', '▓', '█']

        for r in range(HEIGHT):
            for c in range(WIDTH):
                if z_buffer[r, c] != np.inf:
                    depth = z_buffer[r, c]
                    # Normalize depth to [0, 1]
                    normalized_depth = (depth - min_z) / (max_z - min_z) if max_z > min_z else 0.5
                    # Map normalized depth to a shade
                    shade_index = int(normalized_depth * 3.99)
                    screen[r, c] = shades[shade_index]

    # 6. OUTPUT: Print the final image
    print("\nFinal rendered view:")
    for row in screen:
        # Use rstrip() to remove trailing whitespace for cleaner output
        print("".join(row).rstrip())

generate_rotated_torus_view()
<<<C>>>