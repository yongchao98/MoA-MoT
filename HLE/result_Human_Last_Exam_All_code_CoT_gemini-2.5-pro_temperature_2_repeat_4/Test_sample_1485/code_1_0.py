import math
import numpy as np

def solve_torus_rotation():
    """
    This function calculates the new view of a torus after rotation and prints it as ASCII art.
    """
    # Rotation values from the prompt (in degrees)
    rx, ry, rz = 140, 75, 35
    
    # In standard right-handed coordinate systems, positive angles are counter-clockwise.
    # The prompt specifies clockwise rotation, so we use negative angles.
    theta_x = math.radians(-rx)
    theta_y = math.radians(-ry)
    theta_z = math.radians(-rz)

    # Torus and Screen Parameters
    R, r = 2.0, 1.0  # Major and minor radii
    screen_width, screen_height = 50, 22
    u_steps, v_steps = 90, 90 # Resolution of the torus model

    # Screen buffers
    z_buffer = np.full((screen_height, screen_width), float('inf'))
    output_buffer = np.full((screen_height, screen_width), ' ')
    
    # Shading characters: Near -> Light ('░'), Far -> Dark ('█')
    chars = ['░', '▒', '▓', '█']

    # Create rotation matrices
    Rx = np.array([
        [1, 0, 0],
        [0, math.cos(theta_x), -math.sin(theta_x)],
        [0, math.sin(theta_x), math.cos(theta_x)]
    ])
    Ry = np.array([
        [math.cos(theta_y), 0, math.sin(theta_y)],
        [0, 1, 0],
        [-math.sin(theta_y), 0, math.cos(theta_y)]
    ])
    Rz = np.array([
        [math.cos(theta_z), -math.sin(theta_z), 0],
        [math.sin(theta_z), math.cos(theta_z), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation matrix (applied as Rz * Ry * Rx)
    R_combined = Rz @ Ry @ Rx

    # Generate, rotate, and project torus points
    for i in range(u_steps):
        u = 2 * math.pi * i / u_steps
        for j in range(v_steps):
            v = 2 * math.pi * j / v_steps

            # Original point on the torus
            P = np.array([
                (R + r * math.cos(u)) * math.cos(v),
                (R + r * math.cos(u)) * math.sin(v),
                r * math.sin(u)
            ])

            # Rotated point
            P_rotated = R_combined @ P
            x, y, z = P_rotated

            # Projection onto 2D screen
            # Character aspect ratio correction
            y_aspect_ratio_correction = 0.5 
            scale = min(screen_width, screen_height) * 2.5
            screen_x = int(screen_width / 2 + scale * x / (R + r))
            screen_y = int(screen_height / 2 - scale * y_aspect_ratio_correction * y / (R + r))

            # Update Z-buffer if this point is closer to the viewer
            if 0 <= screen_x < screen_width and 0 <= screen_y < screen_height:
                if z < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = z

    # Map z-buffer values to characters for shading
    valid_zs = z_buffer[z_buffer != float('inf')]
    if len(valid_zs) > 0:
        min_z, max_z = valid_zs.min(), valid_zs.max()
        z_range = max_z - min_z
        if z_range == 0: z_range = 1 # Avoid division by zero

        for r_idx in range(screen_height):
            for c_idx in range(screen_width):
                if z_buffer[r_idx, c_idx] != float('inf'):
                    # Normalize depth to range [0, 1]
                    norm_depth = (z_buffer[r_idx, c_idx] - min_z) / z_range
                    # Select character based on normalized depth
                    char_index = min(int(norm_depth * (len(chars))), len(chars) - 1)
                    output_buffer[r_idx, c_idx] = chars[char_index]

    # Print the resulting ASCII art
    for row in output_buffer:
        print("".join(row))
        
    print(f"\nRotation (X, Y, Z) in degrees: ({rx}, {ry}, {rz})")


solve_torus_rotation()