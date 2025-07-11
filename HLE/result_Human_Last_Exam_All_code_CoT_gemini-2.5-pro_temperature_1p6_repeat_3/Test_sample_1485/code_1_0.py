import numpy as np

def generate_rotated_torus_view(rx_deg, ry_deg, rz_deg):
    """
    Generates and prints an ASCII art view of a torus after rotation.
    
    Args:
        rx_deg (float): Rotation in degrees around the X-axis.
        ry_deg (float): Rotation in degrees around the Y-axis.
        rz_deg (float): Rotation in degrees around the Z-axis.
    """
    
    # Helper function to create standard counter-clockwise rotation matrices
    def rotation_matrix(axis, angle_rad):
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        if axis == 'x':
            return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        elif axis == 'y':
            return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        elif axis == 'z':
            return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    # Screen and torus parameters
    screen_width = 48
    screen_height = 22
    R, r = 2.0, 1.0  # Major and minor radii of the torus

    # Convert degrees to radians. Use negative angles for clockwise rotation.
    ax_rad = np.radians(-rx_deg)
    ay_rad = np.radians(-ry_deg)
    az_rad = np.radians(-rz_deg)
    
    # Create individual rotation matrices
    rot_x = rotation_matrix('x', ax_rad)
    rot_y = rotation_matrix('y', ay_rad)
    rot_z = rotation_matrix('z', az_rad)

    # Combine into a single rotation matrix (apply X, then Y, then Z)
    R_total = rot_z @ rot_y @ rot_x

    # Initialize screen grid and z-buffer for depth perception
    screen = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # Generate points on the torus surface
    u_step, v_step = 0.07, 0.07
    for u in np.arange(0, 2 * np.pi, u_step):
        for v in np.arange(0, 2 * np.pi, v_step):
            # Parametric equations for a torus point
            P = np.array([
                (R + r * np.cos(u)) * np.cos(v),
                (R + r * np.cos(u)) * np.sin(v),
                r * np.sin(u)
            ])
            
            # Rotate the point
            P_rot = R_total @ P
            x_rot, y_rot, z_rot = P_rot[0], P_rot[1], P_rot[2]
            
            # Project the 3D point to 2D screen coordinates
            # Scale to fit the screen and adjust for character aspect ratio (approx 2:1 height to width)
            scale = screen_height / 2.8
            screen_x = int(screen_width / 2 + 2 * scale * x_rot)
            screen_y = int(screen_height / 2 - scale * y_rot)
            
            # Z-buffering: if this point is closer than what's already at this pixel, update it
            if 0 <= screen_y < screen_height and 0 <= screen_x < screen_width:
                if z_rot < zbuffer[screen_y][screen_x]:
                    zbuffer[screen_y][screen_x] = z_rot

    # Find the min and max depth values of visible points for shading
    min_z = min(z for row in zbuffer for z in row if z != float('inf'))
    max_z = max(z for row in zbuffer for z in row if z != float('inf'))
    z_range = max_z - min_z
    if z_range == 0: z_range = 1

    # Define the shades from nearest (lightest) to furthest (darkest)
    shades = ['░', '▒', '▓', '█']

    # Fill the screen grid with the appropriate character based on depth
    for y in range(screen_height):
        for x in range(screen_width):
            z = zbuffer[y][x]
            if z != float('inf'):
                # Normalize the depth and select the shade
                normalized_z = (z - min_z) / z_range
                shade_index = min(int(normalized_z * len(shades)), len(shades) - 1)
                screen[y][x] = shades[shade_index]

    # Print the final rendered view
    for row in screen:
        print("".join(row))
        
    # Print the rotation values as requested by the prompt
    print(f"\nRotation applied (X, Y, Z): ({rx_deg}, {ry_deg}, {rz_deg})")

# Run the simulation with the given rotation values
rx, ry, rz = 140, 75, 35
generate_rotated_torus_view(rx, ry, rz)