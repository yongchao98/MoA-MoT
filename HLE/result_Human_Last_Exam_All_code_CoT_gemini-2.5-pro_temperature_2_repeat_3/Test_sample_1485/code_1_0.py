import numpy as np

def render_rotated_torus(rot_x, rot_y, rot_z):
    """
    Generates and prints an ASCII representation of a rotated torus.
    
    The function models a torus, applies a 3D rotation, and renders the result
    using orthographic projection with a z-buffer for depth perception. Shading
    is based on the depth of the surface.
    """

    print(f"Generating view for rotation (X, Y, Z) in degrees: ({rot_x}, {rot_y}, {rot_z})")
    print("-" * 60)

    # --- Parameters ---
    screen_width, screen_height = 48, 22
    R1, R2 = 1.0, 0.5  # Major and minor radii of the torus
    K_scale = screen_width / (2 * (R1 + R2)) * 0.8 # Scale factor for projection
    aspect_ratio_correction = 0.5 # To correct for non-square characters

    # --- Initialization ---
    output_grid = [[' '] * screen_width for _ in range(screen_height)]
    z_buffer = [[float('inf')] * screen_width for _ in range(screen_height)]
    
    # --- Rotation Matrix Setup ---
    # The problem specifies clockwise rotation. Standard matrices are counter-clockwise,
    # so we use negative angles.
    alpha = np.radians(-rot_x)
    beta = np.radians(-rot_y)
    gamma = np.radians(-rot_z)
    
    cos_a, sin_a = np.cos(alpha), np.sin(alpha)
    cos_b, sin_b = np.cos(beta), np.sin(beta)
    cos_g, sin_g = np.cos(gamma), np.sin(gamma)

    # Rotation matrices for each axis
    Rx = np.array([[1, 0, 0], [0, cos_a, -sin_a], [0, sin_a, cos_a]])
    Ry = np.array([[cos_b, 0, sin_b], [0, 1, 0], [-sin_b, 0, cos_b]])
    Rz = np.array([[cos_g, -sin_g, 0], [sin_g, cos_g, 0], [0, 0, 1]])

    # Combined rotation matrix (applied as Rz * Ry * Rx)
    R_combined = Rz @ Ry @ Rx
    
    # --- Generate and Project Points ---
    theta_step = 0.05  # Angle around the major circle
    phi_step = 0.04    # Angle around the minor circle (tube)
    
    for theta in np.arange(0, 2 * np.pi, theta_step):
        for phi in np.arange(0, 2 * np.pi, phi_step):
            # Parametric equations for a standard torus in the xy-plane
            x_std = (R1 + R2 * np.cos(theta)) * np.cos(phi)
            y_std = (R1 + R2 * np.cos(theta)) * np.sin(phi)
            z_std = R2 * np.sin(theta)
            
            # Initial orientation: Rotate torus to be vertical along the Y-axis
            # This is a 90-degree rotation around the X-axis.
            x_init = x_std
            y_init = -z_std
            z_init = y_std
            p_init = np.array([x_init, y_init, z_init])
            
            # Apply the specified rotation
            p_rotated = R_combined @ p_init
            x, y, z = p_rotated[0], p_rotated[1], p_rotated[2]
            
            # Project 3D point to 2D screen
            screen_x = int(screen_width / 2 + K_scale * x)
            screen_y = int(screen_height / 2 - K_scale * y * aspect_ratio_correction)
            
            # Update buffers if the point is visible (closer to the camera)
            if 0 <= screen_y < screen_height and 0 <= screen_x < screen_width:
                if z < z_buffer[screen_y][screen_x]:
                    z_buffer[screen_y][screen_x] = z
    
    # --- Map Depth to Characters for Shading ---
    # Find min and max depth for normalization
    min_z = min(z for row in z_buffer for z in row if z != float('inf'))
    max_z = max(z for row in z_buffer for z in row if z != float('inf'))
    depth_range = max_z - min_z
    
    # Lighter shades (░) are nearer (smaller z)
    shading_chars = ['░', '▒', '▓', '█']

    for r in range(screen_height):
        for c in range(screen_width):
            if z_buffer[r][c] != float('inf'):
                # Normalize depth from 0 to 1
                norm_depth = (z_buffer[r][c] - min_z) / depth_range if depth_range > 0 else 0
                # Select character based on normalized depth
                char_index = int(norm_depth * (len(shading_chars) - 0.001))
                output_grid[r][c] = shading_chars[char_index]
    
    # --- Print Final ASCII Art ---
    for row in output_grid:
        print("".join(row))

if __name__ == '__main__':
    # Rotation values from the problem statement
    rotations_xyz = (140, 75, 35)
    render_rotated_torus(*rotations_xyz)