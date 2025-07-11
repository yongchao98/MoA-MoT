import math

def generate_rotated_torus_view():
    """
    This function generates an ASCII art representation of a torus
    after applying specified 3D rotations. It simulates the process by:
    1. Creating a 3D model of a torus.
    2. Applying clockwise rotations around the X, Y, and Z axes.
    3. Projecting the 3D model onto a 2D screen with depth shading.
    """
    
    # Rotation values from the problem description
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35

    print(f"Calculating view after rotation by ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees around the X, Y, and Z axes...")
    
    # Screen and Torus Parameters
    width, height = 51, 22
    R_major, r_minor = 8.5, 4.5  # Major and minor radii of the torus
    scale = 1.4  # Scale factor for projection
    aspect_ratio = 2.0  # To account for non-square characters in terminals

    # Buffers for rendering
    screen_buffer = [[' '] * width for _ in range(height)]
    z_buffer = [[float('inf')] * width for _ in range(height)]
    
    # Convert degrees to radians. Use negative angles for clockwise rotation.
    ax = -math.radians(rot_x_deg)
    ay = -math.radians(rot_y_deg)
    az = -math.radians(rot_z_deg)
    
    # Pre-calculate sines and cosines for rotation matrices
    cos_ax, sin_ax = math.cos(ax), math.sin(ax)
    cos_ay, sin_ay = math.cos(ay), math.sin(ay)
    cos_az, sin_az = math.cos(az), math.sin(az)

    # Point generation and transformation
    # More steps result in a higher resolution image
    phi_steps = 150
    theta_steps = 80
    
    # Store rotated points to normalize depth later
    rotated_points = []
    
    for i in range(phi_steps):
        phi = 2 * math.pi * i / phi_steps
        cos_phi, sin_phi = math.cos(phi), math.sin(phi)
        
        for j in range(theta_steps):
            theta = 2 * math.pi * j / theta_steps
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            
            # Initial point on the torus (Z-axis is the axis of revolution)
            x = (R_major + r_minor * cos_theta) * cos_phi
            y = (R_major + r_minor * cos_theta) * sin_phi
            z = r_minor * sin_theta
            
            # Apply rotations sequentially: X, then Y, then Z
            # Rotate around X
            y1 = y * cos_ax - z * sin_ax
            z1 = y * sin_ax + z * cos_ax
            x1 = x
            
            # Rotate around Y
            x2 = x1 * cos_ay + z1 * sin_ay
            z2 = -x1 * sin_ay + z1 * cos_ay
            y2 = y1

            # Rotate around Z
            x3 = x2 * cos_az - y2 * sin_az
            y3 = x2 * sin_az + y2 * cos_az
            z3 = z2
            
            rotated_points.append((x3, y3, z3))

    # Determine depth range for shading
    min_z = min(p[2] for p in rotated_points)
    max_z = max(p[2] for p in rotated_points)
    z_range = max_z - min_z
    if z_range == 0: z_range = 1 # Avoid division by zero
    
    shades = ['░', '▒', '▓', '█'] # Near to far

    # Project points onto the 2D screen buffer
    for x_rot, y_rot, z_rot in rotated_points:
        # Orthographic projection with scaling
        sx = int(width / 2 + scale * x_rot)
        sy = int(height / 2 - scale * y_rot / aspect_ratio) # Invert Y for screen coordinates
        
        if 0 <= sy < height and 0 <= sx < width:
            # Z-buffering for occlusion
            if z_rot < z_buffer[sy][sx]:
                z_buffer[sy][sx] = z_rot
                
                # Determine shade based on depth
                shade_index = int(((z_rot - min_z) / z_range) * len(shades))
                shade_index = max(0, min(len(shades) - 1, shade_index))
                screen_buffer[sy][sx] = shades[shade_index]

    # Print the final rendered view
    for row in screen_buffer:
        print("".join(row))

generate_rotated_torus_view()