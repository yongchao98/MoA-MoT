import math

def render_rotated_torus():
    """
    This function generates and prints an ASCII representation of a torus
    after applying specified rotations around the X, Y, and Z axes.
    """
    # The equation to solve: Find the resulting view after the given rotations.
    # The numbers in the equation are the rotation angles.
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35
    
    # --- Configuration ---
    # Screen dimensions
    screen_width = 50
    screen_height = 22
    
    # Torus geometry
    R1 = 1.0  # Minor radius (radius of the tube)
    R2 = 2.0  # Major radius (distance from center to tube's center)
    
    # Viewer and projection settings
    K2 = 5.0  # Distance of the observer from the screen
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) # Scaling factor for projection
    
    # Convert degrees to radians
    rot_x_rad = math.radians(rot_x_deg)
    rot_y_rad = math.radians(rot_y_deg)
    rot_z_rad = math.radians(rot_z_deg)
    
    # Pre-calculate sines and cosines for rotation matrices (clockwise)
    cosA, sinA = math.cos(rot_x_rad), math.sin(rot_x_rad) # For X-axis
    cosB, sinB = math.cos(rot_y_rad), math.sin(rot_y_rad) # For Y-axis
    cosC, sinC = math.cos(rot_z_rad), math.sin(rot_z_rad) # For Z-axis
    
    # Initialize buffers
    output = [[' '] * screen_width for _ in range(screen_height)]
    zbuffer = [[0] * screen_width for _ in range(screen_height)]
    
    # Character set for shading, from light (nearer/brighter) to dark (further/dimmer)
    shading_chars = '░▒▓█'
    
    # --- Main Loop: Iterate through the torus surface ---
    theta = 0
    while theta < 2 * math.pi: # Angle around the tube's cross-section
        phi = 0
        while phi < 2 * math.pi: # Angle around the main hole of the torus
            # --- Point Calculation ---
            costheta, sintheta = math.cos(theta), math.sin(theta)
            cosphi, sinphi = math.cos(phi), math.sin(phi)
            
            # (x, y, z) coords of the point on the torus surface
            circlex = R2 + R1 * costheta
            x = circlex * cosphi
            y = circlex * sinphi
            z = R1 * sintheta
            
            # --- Rotation of the Point (X, then Y, then Z) ---
            # 1. Rotate around X-axis
            x1 = x
            y1 = y * cosA + z * sinA
            z1 = -y * sinA + z * cosA
            
            # 2. Rotate around Y-axis
            x2 = x1 * cosB - z1 * sinB
            y2 = y1
            z2 = x1 * sinB + z1 * cosB
            
            # 3. Rotate around Z-axis
            x_final = x2 * cosC + y2 * sinC
            y_final = -x2 * sinC + y2 * cosC
            z_final = z2
            
            # --- Shading Calculation ---
            # To determine the shade, we calculate the surface normal and its dot product
            # with a light source vector.
            nx = costheta * cosphi
            ny = costheta * sinphi
            nz = sintheta
            
            # Rotate the normal vector just like the point
            nx1, ny1, nz1 = nx, ny * cosA + nz * sinA, -ny * sinA + nz * cosA
            nx2, ny2, nz2 = nx1 * cosB - nz1 * sinB, ny1, nx1 * sinB + nz1 * cosB
            # The final z-component of the normal is nz2
            final_normal_z = nz2
            
            # Luminance is the dot product with light source (0, 0, -1)
            luminance = -final_normal_z
            
            # --- Projection to 2D Screen ---
            z_final_translated = z_final + K2
            if z_final_translated != 0:
                ooz = 1 / z_final_translated # "one over z" for perspective
                
                # Project to screen coordinates
                xp = int(screen_width / 2 + K1 * ooz * x_final)
                yp = int(screen_height / 2 - K1 * ooz * y_final) # Y-axis is inverted
                
                # --- Z-buffering and Output ---
                if 0 <= xp < screen_width and 0 <= yp < screen_height:
                    if ooz > zbuffer[yp][xp]: # If this point is closer than what's stored
                        zbuffer[yp][xp] = ooz
                        
                        # Select character based on luminance
                        if luminance > 0:
                            # Map luminance (0 to 1) to the character set index (0 to 3)
                            char_index = int(luminance * 3.99)
                            char_index = max(0, min(len(shading_chars) - 1, char_index))
                            output[yp][xp] = shading_chars[char_index]
            
            phi += 0.02 # Step through phi
        theta += 0.07 # Step through theta
        
    # --- Final Output ---
    print(f"Applying rotations (X, Y, Z) in degrees: ({rot_x_deg}, {rot_y_deg}, {rot_z_deg})")
    print("Resulting view:")
    print("                                                    ") # Padding to match problem format
    for row in output:
        print("".join(row))

render_rotated_torus()