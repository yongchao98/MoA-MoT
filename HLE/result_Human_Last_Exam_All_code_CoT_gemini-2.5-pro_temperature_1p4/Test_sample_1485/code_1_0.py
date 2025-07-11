import math

def solve():
    """
    This function calculates and prints the 2D ASCII representation of a torus
    after applying a series of 3D rotations.
    """
    
    # Rotation angles in degrees
    x_rot_deg = 140
    y_rot_deg = 75
    z_rot_deg = 35

    print(f"Applying rotation (X, Y, Z): ({x_rot_deg}, {y_rot_deg}, {z_rot_deg}) degrees")
    print("")

    # Convert degrees to radians
    x_rot = math.radians(x_rot_deg)
    y_rot = math.radians(y_rot_deg)
    z_rot = math.radians(z_rot_deg)

    # Torus parameters
    R1 = 2.0  # Major radius
    R2 = 1.0  # Minor radius

    # Screen dimensions
    screen_width = 46
    screen_height = 22

    # Initialize buffers
    # output buffer for the characters
    output = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    # z-buffer for depth checking, initialized to a very large value
    zbuffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # Pre-calculate sines and cosines of rotation angles for efficiency
    cos_x, sin_x = math.cos(x_rot), math.sin(x_rot)
    cos_y, sin_y = math.cos(y_rot), math.sin(y_rot)
    cos_z, sin_z = math.cos(z_rot), math.sin(z_rot)

    # Iterate through the torus surface using two angles, theta and phi
    theta_step = 0.04
    phi_step = 0.02
    
    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            # Parametric equations for a torus with its hole along the Y-axis
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            
            x = (R1 + R2 * cos_theta) * cos_phi
            y = R2 * sin_theta
            z = (R1 + R2 * cos_theta) * sin_phi

            # --- Apply rotations ---
            # 1. Rotate around X-axis
            y_rx = y * cos_x - z * sin_x
            z_rx = y * sin_x + z * cos_x
            x_rx = x
            
            # 2. Rotate around Y-axis
            x_ry = x_rx * cos_y + z_rx * sin_y
            z_ry = -x_rx * sin_y + z_rx * cos_y
            y_ry = y_rx
            
            # 3. Rotate around Z-axis
            x_rz = x_ry * cos_z - y_ry * sin_z
            y_rz = x_ry * sin_z + y_ry * cos_z
            z_rz = z_ry
            
            final_x, final_y, final_z = x_rz, y_rz, z_rz

            # --- Project onto 2D screen ---
            # Orthographic projection
            scale_factor = 6.5
            x_proj = int(screen_width / 2 + scale_factor * final_x)
            # Y is inverted on most screen coordinate systems
            y_proj = int(screen_height / 2 - scale_factor * final_y) 

            # Check if the projected point is on screen and closer than what's already there
            if 0 <= y_proj < screen_height and 0 <= x_proj < screen_width:
                if final_z < zbuffer[y_proj][x_proj]:
                    zbuffer[y_proj][x_proj] = final_z
                    
                    # --- Determine shading based on depth ---
                    # The prompt states: "lighter shades...are nearer"
                    # and "Z-axis...inwards (positive) and outwards (negative)".
                    # So, more negative Z is nearer and should get a lighter shade.
                    shades = '░▒▓█'
                    # Normalize z value to map it to a shade index.
                    # The z-value range for a torus of R1,R2 is approx -(R1+R2) to +(R1+R2)
                    z_range = R1 + R2
                    norm_z = (final_z + z_range) / (2 * z_range)
                    shade_index = int(norm_z * len(shades))
                    # Clamp index to be within bounds of the shades string
                    shade_index = max(0, min(shade_index, len(shades) - 1))
                    
                    output[y_proj][x_proj] = shades[shade_index]

            phi += phi_step
        theta += theta_step

    # Print the final rendered image
    for row in output:
        print("".join(row))

solve()
<<<B>>>