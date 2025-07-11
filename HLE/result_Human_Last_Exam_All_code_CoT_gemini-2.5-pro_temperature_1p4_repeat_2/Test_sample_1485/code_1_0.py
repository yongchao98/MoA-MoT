import math

def solve_torus_rotation():
    """
    Generates and prints an ASCII representation of a torus rotated by given angles.
    This function models a torus, applies 3D rotations, projects it to a 2D screen,
    and calculates shading to produce an ASCII art image.
    """
    # 1. Define rotation angles from the problem statement
    x_rot_deg = 140
    y_rot_deg = 75
    z_rot_deg = 35

    # Convert angles to radians for use in trigonometric functions
    ax = math.radians(x_rot_deg)
    ay = math.radians(y_rot_deg)
    az = math.radians(z_rot_deg)

    # Pre-calculate sines and cosines for efficiency
    sin_ax, cos_ax = math.sin(ax), math.cos(ax)
    sin_ay, cos_ay = math.sin(ay), math.cos(ay)
    sin_az, cos_az = math.sin(az), math.cos(az)

    # 2. Define geometry and screen parameters
    screen_width = 54
    screen_height = 24
    R1 = 1  # Minor radius (thickness of the tube)
    R2 = 2  # Major radius (distance from center of torus to center of tube)
    K2 = 5  # Distance from the viewer to the object
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) # Projection scaling factor

    # Initialize the output screen and a Z-buffer for handling depth
    output = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[-1e9 for _ in range(screen_width)] for _ in range(screen_height)]
    
    # Define shading characters from light to dark
    shading_chars = "░▒▓█"

    # 3. Iterate through points on the torus surface
    theta_step = 0.07
    phi_step = 0.07
    
    theta = 0
    while theta < 2 * math.pi:
        cos_theta, sin_theta = math.cos(theta), math.sin(theta)
        phi = 0
        while phi < 2 * math.pi:
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            
            # Calculate original 3D coordinates of the point
            circlex = R2 + R1 * cos_theta
            x = circlex * cos_phi
            y = circlex * sin_phi
            z = R1 * sin_theta
            
            # 4. Apply rotations sequentially
            # Rotation around X-axis
            y1 = y * cos_ax - z * sin_ax
            z1 = y * sin_ax + z * cos_ax
            x1 = x
            
            # Rotation around Y-axis
            x2 = x1 * cos_ay + z1 * sin_ay
            z2 = -x1 * sin_ay + z1 * cos_ay
            y2 = y1

            # Rotation around Z-axis
            x3 = x2 * cos_az - y2 * sin_az
            y3 = x2 * sin_az + y2 * cos_az
            z3 = z2
            
            # Translate the object away from the camera
            z3 += K2
            
            # 5. Project the 3D point onto the 2D screen
            if z3 == 0: continue
            ooz = 1 / z3
            xp = int(screen_width / 2 + K1 * ooz * x3)
            yp = int(screen_height / 2 - K1 * ooz * y3)

            # 6. Calculate shading based on lighting
            # Rotate the surface normal vector
            nx, ny, nz = cos_theta * cos_phi, cos_theta * sin_phi, sin_theta
            ny1, nz1 = ny * cos_ax - nz * sin_ax, ny * sin_ax + nz * cos_ax
            nx1 = nx
            nx2, nz2 = nx1 * cos_ay + nz1 * sin_ay, -nx1 * sin_ay + nz1 * cos_ay
            ny2 = ny1
            nx3, ny3 = nx2 * cos_az - ny2 * sin_az, nx2 * sin_az + ny2 * cos_az
            nz3 = nz2
            
            # Calculate luminance (dot product with light source at (0,0,-1))
            luminance = -nz3

            # 7. Update buffers if the point is the closest so far
            if 0 <= xp < screen_width and 0 <= yp < screen_height:
                if ooz > zbuffer[yp][xp]:
                    zbuffer[yp][xp] = ooz
                    if luminance > 0:
                        char_index = int(luminance * (len(shading_chars) - 0.01))
                        output[yp][xp] = shading_chars[char_index]
            
            phi += phi_step
        theta += theta_step

    # 8. Print the final result
    print("Executing rotation with the following values:")
    print(f"Rotation around X-axis: {x_rot_deg} degrees")
    print(f"Rotation around Y-axis: {y_rot_deg} degrees")
    print(f"Rotation around Z-axis: {z_rot_deg} degrees")
    print("\nResulting view:\n")
    
    for row in output:
        print("".join(row))

solve_torus_rotation()
<<<B>>>