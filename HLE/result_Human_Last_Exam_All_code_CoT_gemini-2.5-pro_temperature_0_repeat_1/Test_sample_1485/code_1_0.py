import numpy as np
import math

def solve_torus_rotation():
    """
    This function simulates the rotation of a 3D torus and renders the result
    as ASCII art, matching one of the provided answer choices.
    """
    # --- Parameters ---
    # Torus geometry
    R, r = 2.0, 1.0  # Major and minor radii of the torus

    # Screen dimensions for the ASCII output
    screen_width, screen_height = 50, 22

    # Rotation angles in degrees
    theta_x_deg, theta_y_deg, theta_z_deg = 140, 75, 35
    
    print(f"Applying rotation (X, Y, Z) = ({theta_x_deg}, {theta_y_deg}, {theta_z_deg}) degrees")
    print("\n")

    # Convert angles to radians. Use negative angles for clockwise rotation
    # with standard counter-clockwise rotation matrices.
    theta_x = math.radians(-theta_x_deg)
    theta_y = math.radians(-theta_y_deg)
    theta_z = math.radians(-theta_z_deg)

    # --- Rotation Matrices ---
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
    R_total = Rz @ Ry @ Rx

    # --- Buffers ---
    # Initialize the screen with blank spaces
    screen = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    # Initialize the depth buffer with infinity
    z_buffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # --- Simulation Loop ---
    # Iterate through the angles u and v to generate points on the torus surface
    u_step, v_step = 0.05, 0.05
    u = 0
    while u < 2 * math.pi:
        v = 0
        while v < 2 * math.pi:
            # Parametric equations for a point on the torus
            point = np.array([
                (R + r * math.cos(v)) * math.cos(u),
                (R + r * math.cos(v)) * math.sin(u),
                r * math.sin(v)
            ])
            
            # Apply the combined rotation
            rotated_point = R_total @ point
            px, py, pz = rotated_point[0], rotated_point[1], rotated_point[2]

            # Project the 3D point to 2D screen coordinates
            scale = 6
            aspect_ratio_correction = 2.2  # To account for non-square character cells
            screen_x = int(screen_width / 2 + scale * aspect_ratio_correction * px)
            screen_y = int(screen_height / 2 - scale * py) # Y is inverted in screen coords

            # Z-buffer check for occlusion
            if 0 <= screen_x < screen_width and 0 <= screen_y < screen_height:
                if pz < z_buffer[screen_y][screen_x]:
                    z_buffer[screen_y][screen_x] = pz
            v += v_step
        u += u_step

    # --- Shading and Rendering ---
    # Find the min and max depth values for normalization
    min_z, max_z = float('inf'), float('-inf')
    for row in z_buffer:
        for z in row:
            if z != float('inf'):
                min_z = min(min_z, z)
                max_z = max(max_z, z)

    # Define shades from nearest (lightest) to farthest (darkest)
    shades = ['░', '▒', '▓', '█']
    z_range = max_z - min_z
    if z_range == 0: z_range = 1 # Avoid division by zero

    # Fill the screen buffer with the appropriate shade character
    for y in range(screen_height):
        for x in range(screen_width):
            if z_buffer[y][x] != float('inf'):
                # Normalize depth to a 0-1 range
                normalized_z = (z_buffer[y][x] - min_z) / z_range
                # Map the normalized depth to a shade index
                shade_index = int(normalized_z * (len(shades) - 0.001))
                screen[y][x] = shades[shade_index]

    # --- Print the final result ---
    for row in screen:
        print("".join(row).rstrip())

solve_torus_rotation()