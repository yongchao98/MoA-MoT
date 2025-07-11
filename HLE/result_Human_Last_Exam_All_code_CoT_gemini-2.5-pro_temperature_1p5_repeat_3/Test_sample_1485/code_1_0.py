import math

def solve_torus_rotation():
    """
    This function models a torus, rotates it in 3D space,
    and renders a 2D projection using ASCII characters.
    """
    # Rotation angles in degrees
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35

    # Screen and torus parameters
    screen_width = 46
    screen_height = 22
    R1 = 2.0  # Major radius
    R2 = 1.0  # Minor radius
    K2 = 5.0  # Distance of observer from screen
    # K1 is a scaling factor for the projection
    K1 = screen_height * K2 * 3 / (8 * (R1 + R2))

    # Convert degrees to radians
    theta_x = math.radians(rot_x_deg)
    theta_y = math.radians(rot_y_deg)
    theta_z = math.radians(rot_z_deg)

    # Pre-calculate sines and cosines for rotation matrices
    sin_x, cos_x = math.sin(theta_x), math.cos(theta_x)
    sin_y, cos_y = math.sin(theta_y), math.cos(theta_y)
    sin_z, cos_z = math.sin(theta_z), math.cos(theta_z)

    # Initialize buffers for rendering
    # output stores the characters, zbuffer stores the depth of each character
    output = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[0 for _ in range(screen_width)] for _ in range(screen_height)]
    
    shading_chars = "░▒▓█"

    # Iterate through the torus surface using angles theta and phi
    theta = 0
    while theta < 2 * math.pi:
        theta += 0.07  # Step for the major circle
        phi = 0
        while phi < 2 * math.pi:
            phi += 0.02  # Step for the minor circle

            # --- 1. Calculate point and normal vector on original torus ---
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            
            circle_x = R1 + R2 * cos_theta
            
            # Original 3D coordinates of the point
            x = circle_x * cos_phi
            y = circle_x * sin_phi
            z = R2 * sin_theta
            
            # Original normal vector for lighting calculation
            nx = cos_theta * cos_phi
            ny = cos_theta * sin_phi
            nz = sin_theta

            # --- 2. Rotate point and normal vector ---
            # Rotate around X-axis
            y_rot1 = y * cos_x - z * sin_x
            z_rot1 = y * sin_x + z * cos_x
            ny_rot1 = ny * cos_x - nz * sin_x
            nz_rot1 = ny * sin_x + nz * cos_x
            
            # Rotate around Y-axis
            x_rot2 = x * cos_y + z_rot1 * sin_y
            z_rot2 = -x * sin_y + z_rot1 * cos_y
            nx_rot2 = nx * cos_y + nz_rot1 * sin_y
            nz_rot2 = -nx * sin_y + nz_rot1 * cos_y

            # Rotate around Z-axis
            x_rot3 = x_rot2 * cos_z - y_rot1 * sin_z
            y_rot3 = x_rot2 * sin_z + y_rot1 * cos_z
            nx_rot3 = nx_rot2 * cos_z - ny_rot1 * sin_z
            ny_rot3 = nx_rot2 * sin_z + ny_rot1 * cos_z
            
            # Final rotated coordinates and normal vector
            final_x, final_y, final_z = x_rot3, y_rot3, z_rot2
            final_nx, final_ny, final_nz = nx_rot3, ny_rot3, nz_rot2

            # --- 3. Project to 2D screen ---
            z_from_viewer = final_z + K2
            ooz = 1 / z_from_viewer if z_from_viewer != 0 else 0

            # Project to screen coordinates (xp, yp)
            xp = int(screen_width / 2 + K1 * ooz * final_x)
            yp = int(screen_height / 2 - K1 * ooz * final_y)

            # --- 4. Z-buffering and Shading ---
            if 0 <= xp < screen_width and 0 <= yp < screen_height:
                if ooz > zbuffer[yp][xp]:
                    zbuffer[yp][xp] = ooz
                    
                    # Calculate luminance: dot product of normal with light vector
                    # Light source is from above-right-front (0.5, 0.5, -1) for good contrast
                    light_vec = (0.5, 0.5, -1)
                    light_mag = math.sqrt(light_vec[0]**2 + light_vec[1]**2 + light_vec[2]**2)
                    norm_light = tuple(v / light_mag for v in light_vec)
                    
                    luminance = (final_nx * norm_light[0] + 
                                 final_ny * norm_light[1] + 
                                 final_nz * norm_light[2])
                    
                    # Choose character based on luminance
                    if luminance > 0:
                        char_index = int(luminance * len(shading_chars))
                        if char_index >= len(shading_chars):
                            char_index = len(shading_chars) - 1
                        # Lighter shades for higher luminance
                        output[yp][xp] = shading_chars[len(shading_chars) - 1 - char_index]

    # Print the equation and the final rendered view
    print(f"Applying rotations: X={rot_x_deg}, Y={rot_y_deg}, Z={rot_z_deg}")
    print("-" * screen_width)
    for row in output:
        print("".join(row))
    print("-" * screen_width)


# Execute the simulation
solve_torus_rotation()
<<<B>>>