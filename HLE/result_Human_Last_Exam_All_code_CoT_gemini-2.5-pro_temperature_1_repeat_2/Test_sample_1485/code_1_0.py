import math

def solve_torus_rotation():
    """
    This function calculates and prints the 2D ASCII representation of a torus
    after a specified 3D rotation.
    """
    # Input rotation angles in degrees (X, Y, Z)
    x_rot_deg, y_rot_deg, z_rot_deg = 140, 75, 35

    print(f"Calculating view for rotation (X, Y, Z): ({x_rot_deg}, {y_rot_deg}, {z_rot_deg}) degrees")

    # --- 1. Setup Parameters ---

    # Torus geometry
    R1 = 2.0  # Major radius (center of torus to center of tube)
    R2 = 1.0  # Minor radius (radius of the tube)

    # Screen and projection settings
    screen_width = 50
    screen_height = 22
    K2 = 5.0  # Distance from observer to the object
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) # Scaling factor

    # Convert rotation angles to radians.
    # The problem defines positive rotation as clockwise (CW).
    # Standard math libraries use counter-clockwise (CCW) for positive angles.
    # Therefore, we negate the angles to simulate CW rotation.
    A = math.radians(-x_rot_deg)
    B = math.radians(-y_rot_deg)
    C = math.radians(-z_rot_deg)

    # Pre-calculate sines and cosines for rotation matrices
    cosA, sinA = math.cos(A), math.sin(A)
    cosB, sinB = math.cos(B), math.sin(B)
    cosC, sinC = math.cos(C), math.sin(C)

    # Initialize output buffers
    # `output` stores the final characters
    # `z_buffer` stores the inverse depth (1/z) for occlusion
    # `z_values` stores the actual depth for shading
    output = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    z_buffer = [[-float('inf') for _ in range(screen_width)] for _ in range(screen_height)]
    z_values_for_shading = [[None for _ in range(screen_width)] for _ in range(screen_height)]

    # --- 2. Generate, Rotate, and Project Points ---

    # Iterate through the angles of the torus parameterization
    theta_step = 0.07
    phi_step = 0.02
    theta = 0
    while theta < 2 * math.pi:
        costheta, sintheta = math.cos(theta), math.sin(theta)
        phi = 0
        while phi < 2 * math.pi:
            cosphi, sinphi = math.cos(phi), math.sin(phi)

            # (x0, y0, z0): coordinates of a point on the torus surface (initial orientation)
            x0 = (R1 + R2 * cosphi) * costheta
            y0 = (R1 + R2 * cosphi) * sintheta
            z0 = R2 * sinphi

            # Rotate the point sequentially: X -> Y -> Z
            # Rotate around X-axis
            x1 = x0
            y1 = y0 * cosA - z0 * sinA
            z1 = y0 * sinA + z0 * cosA
            # Rotate around Y-axis
            x2 = x1 * cosB + z1 * sinB
            y2 = y1
            z2 = -x1 * sinB + z1 * cosB
            # Rotate around Z-axis
            x_final = x2 * cosC - y2 * sinC
            y_final = x2 * sinC + y2 * cosC
            z_final = z2

            # Project the 3D point to 2D screen coordinates
            z_proj = z_final + K2
            ooz = 1 / z_proj if z_proj != 0 else 0  # One over Z (inverse depth)

            xp = int(screen_width / 2 + K1 * ooz * x_final)
            yp = int(screen_height / 2 - K1 * ooz * y_final)

            # Z-buffer test for occlusion
            if 0 <= xp < screen_width and 0 <= yp < screen_height:
                if ooz > z_buffer[yp][xp]:
                    z_buffer[yp][xp] = ooz
                    z_values_for_shading[yp][xp] = z_final # Store depth for shading
            phi += phi_step
        theta += theta_step

    # --- 3. Determine Shading and Render Final Image ---

    # Find min and max depth of visible points for normalization
    min_z = float('inf')
    max_z = -float('inf')
    for r in range(screen_height):
        for c in range(screen_width):
            if z_values_for_shading[r][c] is not None:
                z_val = z_values_for_shading[r][c]
                min_z = min(min_z, z_val)
                max_z = max(max_z, z_val)
    
    z_range = max_z - min_z if max_z > min_z else 1.0
    
    # Define shading characters from nearest (light) to farthest (dark)
    shade_chars = ['░', '▒', '▓', '█']

    # Assign character based on normalized depth
    for r in range(screen_height):
        for c in range(screen_width):
            if z_values_for_shading[r][c] is not None:
                z_val = z_values_for_shading[r][c]
                # Normalize depth from 0 (nearest) to 1 (farthest)
                norm_z = (z_val - min_z) / z_range
                # Map normalized depth to a character index
                char_index = int(norm_z * (len(shade_chars) - 0.001))
                output[r][c] = shade_chars[char_index]

    # --- 4. Print the Result ---
    print("\n" + "="*screen_width)
    for row in output:
        print("".join(row))
    print("="*screen_width + "\n")

solve_torus_rotation()
<<<B>>>