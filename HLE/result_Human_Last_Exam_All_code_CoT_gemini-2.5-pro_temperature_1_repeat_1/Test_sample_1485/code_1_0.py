import numpy as np

def solve_torus_rotation():
    """
    Calculates and prints the 2D view of a rotated torus.
    """
    # 1. Define rotation parameters
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35
    print(f"Rotation angles in degrees (X, Y, Z): ({rot_x_deg}, {rot_y_deg}, {rot_z_deg})")
    print("-" * 30)

    # Convert degrees to radians for numpy functions
    ax, ay, az = np.radians(rot_x_deg), np.radians(rot_y_deg), np.radians(rot_z_deg)

    # 2. Create clockwise rotation matrices
    # A positive rotation is clockwise. Standard rotation matrices are counter-clockwise.
    # We can use the clockwise formulas directly.
    # Rx_cw(a) = [[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]]
    # Ry_cw(a) = [[cos(a), 0, -sin(a)], [0, 1, 0], [sin(a), 0, cos(a)]]
    # Rz_cw(a) = [[cos(a), sin(a), 0], [-sin(a), cos(a), 0], [0, 0, 1]]
    
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(ax), np.sin(ax)],
        [0, -np.sin(ax), np.cos(ax)]
    ])
    Ry = np.array([
        [np.cos(ay), 0, -np.sin(ay)],
        [0, 1, 0],
        [np.sin(ay), 0, np.cos(ay)]
    ])
    Rz = np.array([
        [np.cos(az), np.sin(az), 0],
        [-np.sin(az), np.cos(az), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation matrix (applied as Rx, then Ry, then Rz)
    rotation_matrix = Rz @ Ry @ Rx

    # 3. Setup rendering canvas
    width, height = 50, 22
    # Torus geometry
    R, r = 2.0, 1.0  # Major and minor radii
    
    # Screen and depth buffers
    screen = np.full((height, width), ' ')
    z_buffer = np.full((height, width), np.inf)
    
    # Shading based on depth (near to far)
    chars = "░▒▓█"
    
    # Projection scaling and offset
    # Character aspect ratio (height/width) is roughly 2, so we scale y by 0.5
    x_scale = width / (2 * (R + r)) * 0.9
    y_scale = height / (2 * (R + r)) * 0.9 * 0.5 
    x_offset = width / 2
    y_offset = height / 2

    # 4. Generate, rotate, and project points
    theta_step = 0.04  # Steps around the major circle
    phi_step = 0.04    # Steps around the minor circle (tube)

    for theta in np.arange(0, 2 * np.pi, theta_step):
        for phi in np.arange(0, 2 * np.pi, phi_step):
            # Parametric equations for a torus with its hole along the Y-axis
            x0 = (R + r * np.cos(phi)) * np.cos(theta)
            y0 = r * np.sin(phi)
            z0 = (R + r * np.cos(phi)) * np.sin(theta)
            
            p0 = np.array([x0, y0, z0])
            p_rot = rotation_matrix @ p0
            x, y, z = p_rot[0], p_rot[1], p_rot[2]

            # Project to 2D screen coordinates
            # Y-axis runs bottom-to-top, so we invert y for screen's top-to-bottom indexing
            screen_x = int(x * x_scale + x_offset)
            screen_y = int(-y * y_scale + y_offset)

            # Z-buffer check for occlusion
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if z < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = z

    # 5. Determine shading and fill screen buffer
    visible_z = z_buffer[z_buffer != np.inf]
    if len(visible_z) > 0:
        min_z, max_z = np.min(visible_z), np.max(visible_z)
        z_range = max_z - min_z if max_z > min_z else 1.0

        for y_idx in range(height):
            for x_idx in range(width):
                if z_buffer[y_idx, x_idx] != np.inf:
                    # Normalize depth to [0, 1]
                    norm_z = (z_buffer[y_idx, x_idx] - min_z) / z_range
                    # Select character based on normalized depth
                    char_index = min(int(norm_z * (len(chars))), len(chars) - 1)
                    screen[y_idx, x_idx] = chars[char_index]

    # Print the final rendered view
    for row in screen:
        print("".join(row))

solve_torus_rotation()