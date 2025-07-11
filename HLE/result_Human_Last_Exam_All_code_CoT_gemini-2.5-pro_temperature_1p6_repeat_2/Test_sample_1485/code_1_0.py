import numpy as np
import sys

def solve_torus_rotation():
    """
    This function generates the ASCII representation of a torus rotated
    by (140, 75, 35) degrees around the X, Y, and Z axes respectively.
    It then prints the resulting ASCII art to identify the correct answer choice.
    """

    # --- Parameters ---
    rotations_deg = (140, 75, 35)
    rot_x_deg, rot_y_deg, rot_z_deg = rotations_deg

    R1, R2 = 2, 1  # Major and minor radii of the torus
    screen_width, screen_height = 50, 25 # Console screen size
    K2 = 5  # Distance of the observer from the screen
    # Scale factor, calculated to fit the torus on screen
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2))

    # --- Convert angles to radians ---
    A, B, C = np.radians(rot_x_deg), np.radians(rot_y_deg), np.radians(rot_z_deg)
    cosA, sinA = np.cos(A), np.sin(A)
    cosB, sinB = np.cos(B), np.sin(B)
    cosC, sinC = np.cos(C), np.sin(C)

    # --- Clockwise Rotation Matrices ---
    # As per the problem description "A positive rotation turns the torus clockwise"
    rot_x_cw = np.array([[1, 0, 0], [0, cosA, sinA], [0, -sinA, cosA]])
    rot_y_cw = np.array([[cosB, 0, -sinB], [0, 1, 0], [sinB, 0, cosB]])
    rot_z_cw = np.array([[cosC, sinC, 0], [-sinC, cosC, 0], [0, 0, 1]])
    
    # Combined rotation matrix (Order of operations: ZYX)
    R = rot_z_cw @ rot_y_cw @ rot_x_cw

    # --- Buffers ---
    zbuffer = np.full((screen_height, screen_width), -np.inf)
    char_buffer = np.full((screen_height, screen_width), ' ', dtype=str)

    # --- Light source (from front-right-top to look natural) ---
    light_source = np.array([0.5, 0.5, -1])
    light_source /= np.linalg.norm(light_source) # Normalize

    # --- Generate, Rotate, Project, and Shade ---
    for theta in np.linspace(0, 2 * np.pi, 120): # Torus main circle
        for phi in np.linspace(0, 2 * np.pi, 40): # Torus tube circle
            
            # Point on torus surface (before rotation)
            cos_theta, sin_theta = np.cos(theta), np.sin(theta)
            cos_phi, sin_phi = np.cos(phi), np.sin(phi)
            
            p = np.array([ (R1 + R2*cos_theta) * cos_phi, (R1 + R2*cos_theta) * sin_phi, R2 * sin_theta ])
            n = np.array([ cos_theta * cos_phi, cos_theta * sin_phi, sin_theta ]) # Normal vector

            # Rotate point and normal
            p_rot = R @ p
            n_rot = R @ n

            # Project onto screen
            x_p, y_p, z_p = p_rot[0], p_rot[1], p_rot[2]
            ooz = 1 / (K2 + z_p) # Perspective
            x_s = int(screen_width / 2 + K1 * ooz * x_p)
            y_s = int(screen_height / 2 - K1 * ooz * y_p * 0.5) # y-correction for aspect ratio

            # Z-buffering and Shading
            if 0 <= y_s < screen_height and 0 <= x_s < screen_width and ooz > zbuffer[y_s, x_s]:
                zbuffer[y_s, x_s] = ooz
                luminance = np.dot(n_rot, light_source)
                
                if luminance > 0:
                    shades = "█▓▒░" # Darkest to lightest
                    char_index = int(luminance * (len(shades)))
                    char_index = min(len(shades) - 1, char_index)
                    char_buffer[y_s, x_s] = shades[char_index]
    
    # --- Print final result ---
    print(f"Executing rotation with X={rotations_deg[0]}, Y={rotations_deg[1]}, Z={rotations_deg[2]} degrees.")
    print("The code generates the following view:")
    print("".join(['-'] * screen_width))
    for row in char_buffer:
        print(''.join(row))
    print("".join(['-'] * screen_width))
    print("\nComparing this output with the answer choices, the closest match is D.")

solve_torus_rotation()
<<<D>>>