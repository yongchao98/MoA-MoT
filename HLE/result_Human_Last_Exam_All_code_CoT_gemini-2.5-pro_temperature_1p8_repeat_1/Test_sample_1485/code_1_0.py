import math

def generate_torus_view():
    """
    This function generates an ASCII representation of a torus after applying
    specified rotations around the X, Y, and Z axes.
    """
    # --- Configuration ---
    width, height = 52, 22
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35
    
    # Use negative angles for clockwise rotation with standard counter-clockwise matrices
    angle_x = -rot_x_deg * math.pi / 180
    angle_y = -rot_y_deg * math.pi / 180
    angle_z = -rot_z_deg * math.pi / 180

    # Torus geometry
    R1, R2 = 15, 6  # Major and minor radii

    # Buffers for rendering
    zbuffer = [[float('inf')] * width for _ in range(height)]
    output = [[' '] * width for _ in range(height)]

    # Pre-calculate trigonometric values for rotation matrices
    cos_ax, sin_ax = math.cos(angle_x), math.sin(angle_x)
    cos_ay, sin_ay = math.cos(angle_y), math.sin(angle_y)
    cos_az, sin_az = math.cos(angle_z), math.sin(angle_z)
    
    # --- Step 1 & 2: Generate points on the torus and rotate them ---
    # Iterate over the two angles that define the torus surface
    theta = 0
    while theta < 2 * math.pi:
        theta += 0.07  # Step size for the major circle
        cos_theta, sin_theta = math.cos(theta), math.sin(theta)
        
        phi = 0
        while phi < 2 * math.pi:
            phi += 0.04  # Step size for the minor circle
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            
            # Parametric equations for a torus in the xy-plane
            circlex = R1 + R2 * cos_theta
            x0 = circlex * cos_phi
            y0 = circlex * sin_phi
            z0 = R2 * sin_theta

            # Rotate point around the X-axis
            y1 = y0 * cos_ax - z0 * sin_ax
            z1 = y0 * sin_ax + z0 * cos_ax
            x1 = x0
            
            # Rotate point around the Y-axis
            x2 = x1 * cos_ay + z1 * sin_ay
            z2 = -x1 * sin_ay + z1 * cos_ay
            y2 = y1

            # Rotate point around the Z-axis
            x3 = x2 * cos_az - y2 * sin_az
            y3 = x2 * sin_az + y2 * cos_az
            z3 = z2
            
            # --- Step 3 & 4: Project to 2D and use Z-buffer ---
            # Orthographic projection onto the screen
            x_proj = int(width / 2 + x3)
            y_proj = int(height / 2 - y3) # Y is inverted on screen coordinates

            # Z-buffering check: only draw if the point is closer to the screen
            if 0 <= x_proj < width and 0 <= y_proj < height:
                if z3 < zbuffer[y_proj][x_proj]:
                    zbuffer[y_proj][x_proj] = z3
    
    # --- Step 5: Shading based on depth ---
    # Find the min and max depth of the visible points
    min_z = float('inf')
    max_z = float('-inf')
    for r in range(height):
        for c in range(width):
            if zbuffer[r][c] != float('inf'):
                min_z = min(min_z, zbuffer[r][c])
                max_z = max(max_z, zbuffer[r][c])

    # Lighter shades ('░') for nearer points (smaller z), darker ('█') for farther
    shades = "░▒▓█"
    if max_z > min_z:
        z_range = max_z - min_z
        for r in range(height):
            for c in range(width):
                if zbuffer[r][c] != float('inf'):
                    # Normalize z-value to [0, 1] range
                    norm_z = (zbuffer[r][c] - min_z) / z_range
                    # Map normalized value to a shade index
                    shade_index = int(norm_z * (len(shades) - 0.001)) # -0.001 to prevent index out of bounds
                    output[r][c] = shades[shade_index]
    
    # --- Step 6: Print the final result ---
    print(f"Result for rotation (X={rot_x_deg}°, Y={rot_y_deg}°, Z={rot_z_deg}°):")
    print()
    for row in output:
        print(''.join(row))

if __name__ == '__main__':
    generate_torus_view()
