import numpy as np

def render_torus(rot_x, rot_y, rot_z):
    """
    Renders a rotated torus using ASCII characters.
    """
    # Screen and model parameters
    screen_width = 45
    screen_height = 22
    output = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]
    
    R = 2.0  # Major radius of the torus
    r = 1.0  # Minor radius of the torus
    
    # Character map for shading from brightest to darkest
    shades = '░▒▓█'
    
    # Convert degrees to radians for numpy functions
    angle_x = np.radians(rot_x)
    angle_y = np.radians(rot_y)
    angle_z = np.radians(rot_z)
    
    # Create standard counter-clockwise rotation matrices
    cx, sx = np.cos(angle_x), np.sin(angle_x)
    cy, sy = np.cos(angle_y), np.sin(angle_y)
    cz, sz = np.cos(angle_z), np.sin(angle_z)
    
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    
    # Combine rotations. The order is extrinsic Z, then Y, then X.
    R_total = np.dot(Rz, np.dot(Ry, Rx))
    
    # Light source is from the observer's direction (from negative Z)
    light_vector = np.array([0, 0, -1])

    # Iterate through the surface of the torus
    theta_step = 0.07  # Angle around the major circle
    phi_step = 0.07    # Angle around the minor circle (tube)
    
    theta = 0
    while theta < 2 * np.pi:
        phi = 0
        while phi < 2 * np.pi:
            # Parametric equations for a torus with its hole along the X-axis
            x_initial = r * np.cos(phi)
            y_initial = (R + r * np.sin(phi)) * np.cos(theta)
            z_initial = (R + r * np.sin(phi)) * np.sin(theta)
            p_initial = np.array([x_initial, y_initial, z_initial])

            # Normal vector for the surface point (points outward from the tube)
            nx_initial = np.cos(phi)
            ny_initial = np.sin(phi) * np.cos(theta)
            nz_initial = np.sin(phi) * np.sin(theta)
            n_initial = np.array([nx_initial, ny_initial, nz_initial])

            # Apply the combined rotation to the point and the normal vector
            p_rotated = np.dot(R_total, p_initial)
            n_rotated = np.dot(R_total, n_initial)

            # Calculate luminance based on the light source
            luminance = np.dot(n_rotated, light_vector)

            # Only draw surfaces facing the camera
            if luminance > 0:
                # Project the 3D point to 2D screen coordinates
                # Scale factors are chosen to fit the aspect ratio of terminal characters
                screen_x = int(screen_width / 2 + 2.5 * p_rotated[0])
                screen_y = int(screen_height / 2 - 1.2 * p_rotated[1])

                if 0 <= screen_x < screen_width and 0 <= screen_y < screen_height:
                    # Z-buffer check for occlusion
                    if p_rotated[2] < zbuffer[screen_y][screen_x]:
                        zbuffer[screen_y][screen_x] = p_rotated[2]
                        # Map luminance to a shade character
                        # High luminance -> bright -> low index character
                        shade_index = int(luminance * len(shades))
                        shade_index = min(shade_index, len(shades) - 1)
                        # High luminance = darker char, let's reverse
                        # To match light char=bright
                        shades_reversed = shades
                        output[screen_y][screen_x] = shades_reversed[shade_index]
            phi += phi_step
        theta += theta_step
    
    # Print the rotation values and the final rendered torus
    print(f"Rotation(X={rot_x}, Y={rot_y}, Z={rot_z}):\n")
    for row in output:
        print("".join(row))

# The rotation values from the problem statement
x_rot, y_rot, z_rot = 140, 75, 35
render_torus(x_rot, y_rot, z_rot)