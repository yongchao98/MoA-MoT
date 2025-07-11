import numpy as np

def render_torus(rot_x_deg, rot_y_deg, rot_z_deg):
    """
    Generates and prints an ASCII representation of a rotated torus.
    
    The rendering uses a right-handed coordinate system and counter-clockwise
    rotations, which is standard. To match the problem's left-handed,
    clockwise rotation system, the signs of the angles are inverted.
    """
    # Torus and Screen Parameters
    R1 = 1.0  # Minor radius (tube thickness)
    R2 = 2.0  # Major radius (distance from center of torus to center of tube)
    screen_width = 60
    screen_height = 24
    
    # Observer/camera distance from the screen
    K2 = 5
    # Scale factor for projection
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2))

    # Convert degrees to radians (inverting angles for LH/CW to RH/CCW conversion)
    rot_x = np.radians(-rot_x_deg)
    rot_y = np.radians(-rot_y_deg)
    rot_z = np.radians(-rot_z_deg)

    # Rotation matrices for counter-clockwise rotation in a right-handed system
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rot_x), -np.sin(rot_x)],
        [0, np.sin(rot_x), np.cos(rot_x)]
    ])
    Ry = np.array([
        [np.cos(rot_y), 0, np.sin(rot_y)],
        [0, 1, 0],
        [-np.sin(rot_y), 0, np.cos(rot_y)]
    ])
    Rz = np.array([
        [np.cos(rot_z), -np.sin(rot_z), 0],
        [np.sin(rot_z), np.cos(rot_z), 0],
        [0, 0, 1]
    ])
    # Combined rotation matrix (order: X, then Y, then Z)
    R = Rz @ Ry @ Rx

    # Initialize buffers: one for characters, one for depth (z-buffer)
    output_buffer = np.full((screen_height, screen_width), ' ')
    z_buffer = np.full((screen_height, screen_width), -np.inf)

    # Light source vector (from the observer's direction)
    light_source = np.array([0, 0, -1])

    # Define the density of points to draw on the torus surface
    theta_spacing = 0.07  # For the tube's cross-section
    phi_spacing = 0.02    # For the main ring of the torus

    # Iterate through the torus's surface angles
    for theta in np.arange(0, 2 * np.pi, theta_spacing):
        for phi in np.arange(0, 2 * np.pi, phi_spacing):
            # Calculate the point on the original unrotated torus
            cos_theta, sin_theta = np.cos(theta), np.sin(theta)
            cos_phi, sin_phi = np.cos(phi), np.sin(phi)
            
            x = (R2 + R1 * cos_theta) * cos_phi
            y = (R2 + R1 * cos_theta) * sin_phi
            z = R1 * sin_theta
            point = np.array([x, y, z])

            # Calculate the normal vector at that point
            nx = cos_theta * cos_phi
            ny = cos_theta * sin_phi
            nz = sin_theta
            normal = np.array([nx, ny, nz])

            # Apply the rotation to both the point and its normal vector
            rotated_point = R @ point
            rotated_normal = R @ normal
            
            # Calculate luminance based on the angle between the surface normal and light source
            luminance = np.dot(rotated_normal, light_source)

            # Only draw surfaces that are facing the light source (and the camera)
            if luminance > 0:
                # Project the 3D point to the 2D screen
                x_rot, y_rot, z_rot = rotated_point
                
                # 'ooz' (one over z) is used for perspective projection and z-buffering
                ooz = 1 / (K2 - z_rot)
                
                # Calculate screen coordinates (xp, yp)
                xp = int(screen_width / 2 + K1 * ooz * x_rot)
                # Correct for non-square character aspect ratio (terminals are ~2:1 height:width)
                yp = int(screen_height / 2 + K1 * ooz * y_rot * 0.5)

                # If the point is within the screen bounds and is closer than anything drawn before
                if 0 <= xp < screen_width and 0 <= yp < screen_height and ooz > z_buffer[yp, xp]:
                    # Update the z-buffer with the new closest depth
                    z_buffer[yp, xp] = ooz
                    
                    # Map luminance to a character shade. Brighter surfaces get lighter characters.
                    if luminance > 0.8:
                        char = '░'
                    elif luminance > 0.5:
                        char = '▒'
                    elif luminance > 0.2:
                        char = '▓'
                    else:
                        char = '█'
                    output_buffer[yp, xp] = char
    
    # Print the final rendered image
    print(f"View after rotation by ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees around X, Y, Z axes:\n")
    for row in output_buffer:
        print("".join(row).center(70))


# The rotation angles from the problem
x_rot, y_rot, z_rot = 140, 75, 35
render_torus(x_rot, y_rot, z_rot)