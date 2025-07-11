import numpy as np
import math

def solve_torus_rotation():
    """
    This function generates an ASCII representation of a torus after applying
    a series of rotations.
    """
    
    # Torus parameters
    R, r = 2.0, 1.0  # Major and minor radii

    # Rotation angles in degrees, as specified in the problem
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35

    # Screen display parameters
    screen_width = 50
    screen_height = 25
    canvas = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    z_buffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # Orthographic projection and centering parameters
    scale = 11
    x_shift = screen_width / 2.2
    y_shift = screen_height / 2

    # --- Rotation Matrices (for standard counter-clockwise rotation) ---
    def rotation_matrix(axis, angle_rad):
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        if axis == 'x':
            return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        elif axis == 'y':
            return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        elif axis == 'z':
            return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        return np.identity(3)

    # --- Set up rotations ---
    # 1. Initial orientation: Torus in XZ-plane (a 90-degree rotation around X)
    initial_rot = rotation_matrix('x', math.radians(90))

    # 2. Specified rotations: Positive rotation is clockwise, so use negative angles.
    rot_x = rotation_matrix('x', -math.radians(rot_x_deg))
    rot_y = rotation_matrix('y', -math.radians(rot_y_deg))
    rot_z = rotation_matrix('z', -math.radians(rot_z_deg))

    # 3. Combine rotations: initial, then X, Y, Z
    final_rotation = rot_z @ rot_y @ rot_x @ initial_rot

    # --- Lighting ---
    # Light vector (simulates light from the top-right-front)
    light_vec = np.array([-0.5, 0.5, -0.707])
    light_vec /= np.linalg.norm(light_vec) # Normalize the vector

    # --- Main Loop ---
    # Iterate through the angles of the torus surface
    theta_step, phi_step = 0.05, 0.05
    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            # Parametric equations for a standard torus
            c_th, s_th = math.cos(theta), math.sin(theta)
            c_ph, s_ph = math.cos(phi), math.sin(phi)
            
            x = (R + r * c_ph) * c_th
            y = (R + r * c_ph) * s_th
            z = r * s_ph
            point = np.array([x, y, z])
            
            normal = np.array([c_ph * c_th, c_ph * s_th, s_ph])
            
            # Apply the final combined rotation to the point and its normal
            rotated_point = final_rotation @ point
            rotated_normal = final_rotation @ normal
            
            # Get projected coordinates
            px, py, pz = rotated_point[0], rotated_point[1], rotated_point[2]
            
            # Calculate luminance based on the light source
            luminance = np.dot(rotated_normal, light_vec)
            
            # Render only surfaces facing the light
            if luminance > 0:
                screen_x = int(x_shift + scale * px)
                screen_y = int(y_shift - scale * py) # Screen y is inverted from model y
                
                # Check bounds and Z-buffer
                if 0 <= screen_x < screen_width and 0 <= screen_y < screen_height:
                    if pz < z_buffer[screen_y][screen_x]:
                        z_buffer[screen_y][screen_x] = pz
                        
                        # Map luminance to shade character
                        if luminance > 0.85: char = '░'
                        elif luminance > 0.55: char = '▒'
                        elif luminance > 0.2: char = '▓'
                        else: char = '█'
                        canvas[screen_y][screen_x] = char
            phi += phi_step
        theta += theta_step

    # --- Print Final Output ---
    print(f"Applying rotation (X, Y, Z): ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees")
    print("Final rendered view:")
    
    # Format the output to match the style of the answer choices
    for _ in range(2): print() # Add leading newlines
    for row in canvas:
        print("                       " + "".join(row))
    for _ in range(6): print() # Add trailing newlines

# Run the simulation
solve_torus_rotation()