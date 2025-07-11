import math

def solve():
    """
    This function models a torus, applies 3D rotations, and renders the result
    as ASCII art, according to the problem description.
    """

    # 1. Define rotation parameters as specified in the problem
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35
    
    # As per instructions, output the numbers from the rotation
    print(f"Applying clockwise rotation (X, Y, Z) = ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees.")
    print("")

    # Standard math functions use counter-clockwise rotation, so we negate the angles
    # for the specified clockwise rotation.
    ax = math.radians(-rot_x_deg)
    ay = math.radians(-rot_y_deg)
    az = math.radians(-rot_z_deg)

    # Pre-calculate sines and cosines for efficiency
    ca, sa = math.cos(ax), math.sin(ax)
    cb, sb = math.cos(ay), math.sin(ay)
    cc, sc = math.cos(az), math.sin(az)

    # 2. Setup rendering parameters
    # Screen and torus dimensions
    width, height = 50, 22
    R1, R2 = 1, 2  # Minor and major radii
    x_scale, y_scale = 10, 5 # Scaling factors for screen projection

    # Initialize buffers for rendering
    # Z-buffer stores depth to handle occlusions
    zbuffer = [[-float('inf')] * width for _ in range(height)]
    # Output buffer stores the final characters
    output_buffer = [[' '] * width for _ in range(height)]
    
    # Define a light source vector (normalized) for shading
    light_vec = (0, 1, -1)
    norm_l = math.sqrt(sum(i**2 for i in light_vec))
    lx, ly, lz = (i / norm_l for i in light_vec)

    # 3. Generate, rotate, and project torus points
    u_step, v_step = 0.07, 0.03 # Step size for iterating over torus surface
    u = 0
    while u < 2 * math.pi:
        u += u_step
        cu, su = math.cos(u), math.sin(u)
        v = 0
        while v < 2 * math.pi:
            v += v_step
            cv, sv = math.cos(v), math.sin(v)

            # Point on the original torus surface (oriented upright)
            x0 = (R2 + R1 * cv) * cu
            y0 = R1 * sv
            z0 = (R2 + R1 * cv) * su
            
            # Surface normal on the original torus
            nx0 = cv * cu
            ny0 = sv
            nz0 = cv * su

            # Rotate the point and the normal vector sequentially (Z, then Y, then X)
            # Z-axis rotation
            x_r1, y_r1 = x0 * cc - y0 * sc, x0 * sc + y0 * cc
            nx_r1, ny_r1 = nx0 * cc - ny0 * sc, nx0 * sc + ny0 * cc
            # Y-axis rotation
            x_r2, z_r2 = x_r1 * cb + z0 * sb, -x_r1 * sb + z0 * cb
            nx_r2, nz_r2 = nx_r1 * cb + nz0 * sb, -nx_r1 * sb + nz0 * cb
            # X-axis rotation
            y_r3, z = y_r1 * ca - z_r2 * sa, y_r1 * sa + z_r2 * ca
            ny, nz = ny_r1 * ca - nz_r2 * sa, ny_r1 * sa + nz_r2 * ca
            x, y = x_r2, y_r3
            nx = nx_r2
            
            # Calculate luminance (brightness) using the dot product with the light vector
            luminance = nx * lx + ny * ly + nz * lz

            # Process only points facing the light
            if luminance > 0:
                # Project the 3D point to the 2D screen
                screen_x = int(width / 2 + x_scale * x)
                screen_y = int(height / 2 - y_scale * y)
                
                # Z-buffering check
                if 0 <= screen_x < width and 0 <= screen_y < height:
                    if z > zbuffer[screen_y][screen_x]:
                        zbuffer[screen_y][screen_x] = z
                        # Select character based on luminance
                        chars = "░▒▓█"
                        char_idx = int(luminance * 8)
                        if char_idx > 3: char_idx = 3 # Clamp to max index
                        
                        # Adjusting thresholds to better match the style of the answer choices
                        if luminance > 0.8: char_idx = 3
                        elif luminance > 0.6: char_idx = 2
                        elif luminance > 0.3: char_idx = 1
                        else: char_idx = 0
                        output_buffer[screen_y][screen_x] = chars[char_idx]

    # 4. Print the final rendered image
    for row in output_buffer:
        print("".join(row))

solve()