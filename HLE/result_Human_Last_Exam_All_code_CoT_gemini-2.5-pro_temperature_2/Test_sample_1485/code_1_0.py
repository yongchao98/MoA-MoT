import numpy as np

def solve_torus_rotation():
    """
    This function generates the ASCII view of a torus after specified 3D rotations.
    """
    
    # --- Configuration ---
    # Screen dimensions for the ASCII output
    width, height = 50, 22

    # Torus parameters: R is the major radius (center of torus to center of tube)
    # r is the minor radius (radius of the tube)
    R, r = 2.0, 1.0

    # Rotation angles in degrees. Positive is clockwise as per the problem description.
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35
    
    # Output the equation parameters as requested
    print(f"Applying rotations (X, Y, Z): ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees")
    print("-" * (width + 12))

    # The characters used for shading, from nearest ('░') to farthest ('█')
    shades = "░▒▓█"

    # --- Mathematical Setup ---
    # Convert degrees to radians. Standard rotation matrices are counter-clockwise,
    # so to perform a clockwise rotation, we use the negative of the angle.
    rot_x_rad = np.deg2rad(-rot_x_deg)
    rot_y_rad = np.deg2rad(-rot_y_deg)
    rot_z_rad = np.deg2rad(-rot_z_deg)

    # Create rotation matrices for each axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rot_x_rad), -np.sin(rot_x_rad)],
        [0, np.sin(rot_x_rad),  np.cos(rot_x_rad)]
    ])

    Ry = np.array([
        [np.cos(rot_y_rad), 0, np.sin(rot_y_rad)],
        [0, 1, 0],
        [-np.sin(rot_y_rad), 0, np.cos(rot_y_rad)]
    ])

    Rz = np.array([
        [np.cos(rot_z_rad), -np.sin(rot_z_rad), 0],
        [np.sin(rot_z_rad),  np.cos(rot_z_rad), 0],
        [0, 0, 1]
    ])

    # Combine rotation matrices. The order of application is Z, then Y, then X.
    R_total = Rz @ Ry @ Rx

    # --- Rendering Setup ---
    # Initialize a z-buffer with infinity and a screen buffer with blank spaces
    z_buffer = np.full((height, width), np.inf)
    screen = np.full((height, width), ' ')

    # --- Point Generation and Projection ---
    # Iterate over the torus surface angles (theta and phi) to generate a dense point cloud.
    theta_steps, phi_steps = 150, 150
    for theta in np.linspace(0, 2 * np.pi, theta_steps):
        for phi in np.linspace(0, 2 * np.pi, phi_steps):
            # Parametric equations for a point on the torus surface
            x = (R + r * np.cos(phi)) * np.cos(theta)
            y = (R + r * np.cos(phi)) * np.sin(theta)
            z = r * np.sin(phi)
            
            # Rotate the point using the combined rotation matrix
            rotated_point = R_total @ np.array([x, y, z])
            rx, ry, rz = rotated_point
            
            # Project the 3D point to 2D screen coordinates.
            # Scale controls the size of the torus on screen.
            # Y-axis is scaled by 0.5 to account for character aspect ratio (characters are taller than they are wide).
            scale = 6.0
            screen_x = int(width / 2 + scale * rx)
            screen_y = int(height / 2 - scale * ry * 0.5)

            # Z-buffering: only draw the point if it's within bounds and closer to the observer
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if rz < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = rz

    # --- Shading and Final Rendering ---
    # Find the min and max depth values of the visible points to normalize the depth
    min_z = np.min(z_buffer[z_buffer != np.inf])
    max_z = np.max(z_buffer)

    # For each pixel, determine the correct character based on its depth
    for y in range(height):
        for x in range(width):
            if z_buffer[y, x] != np.inf:
                # Normalize depth from 0.0 (nearest) to 1.0 (farthest)
                depth = (z_buffer[y, x] - min_z) / (max_z - min_z)
                # Map the normalized depth to a character index
                shade_index = int(depth * (len(shades) - 0.001))
                screen[y, x] = shades[shade_index]

    # --- Print the Final ASCII Art ---
    for row in screen:
        print("".join(row))

if __name__ == '__main__':
    solve_torus_rotation()