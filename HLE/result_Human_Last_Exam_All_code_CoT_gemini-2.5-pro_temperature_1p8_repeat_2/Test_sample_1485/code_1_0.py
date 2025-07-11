import numpy as np
import math

def solve_torus_rotation():
    """
    This function simulates the rotation of a 3D torus and renders it as ASCII art.
    """
    # Rotation values from the problem description
    rx, ry, rz = 140, 75, 35

    # --- Step 1: Define rotation matrices ---
    alpha = math.radians(rx)
    beta = math.radians(ry)
    gamma = math.radians(rz)

    Rx = np.array([
        [1, 0, 0],
        [0, math.cos(alpha), -math.sin(alpha)],
        [0, math.sin(alpha), math.cos(alpha)]
    ])
    Ry = np.array([
        [math.cos(beta), 0, math.sin(beta)],
        [0, 1, 0],
        [-math.sin(beta), 0, math.cos(beta)]
    ])
    Rz = np.array([
        [math.cos(gamma), -math.sin(gamma), 0],
        [math.sin(gamma), math.cos(gamma), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation matrix (applying X, then Y, then Z rotation)
    R = Rz @ Ry @ Rx

    # --- Step 2: Setup rendering parameters ---
    screen_width = 50
    screen_height = 25
    z_buffer = np.full((screen_height, screen_width), float('inf'))
    output = np.full((screen_height, screen_width), ' ')

    # Torus geometry
    R1 = 2.0  # Major radius
    R2 = 1.0  # Minor radius
    
    # Projection scaling factor
    screen_scale = 7.0

    # Shading characters: The prompt describes lighter shades as being "nearer".
    # This is achieved by lighting the model from the front (viewer's position).
    # A high luminance means the surface faces the viewer, appearing bright and "near".
    # Characters are from light '░' to dark '█'.
    shades = '░▒▓█'
    light_source = np.array([0, 0, 1]) # From the viewer's direction

    # --- Step 3: Iterate through torus surface and render ---
    theta_step = 0.07  # Around the major circle
    phi_step = 0.035 # Around the minor circle

    for theta in np.arange(0, 2 * math.pi, theta_step):
        for phi in np.arange(0, 2 * math.pi, phi_step):
            # Original coordinates of a point on the torus (in XZ plane)
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)

            P_orig = np.array([
                (R1 + R2 * cos_phi) * cos_theta,
                R2 * sin_phi,
                (R1 + R2 * cos_phi) * sin_theta
            ])
            
            # Original surface normal vector (pointing outwards)
            N_orig = np.array([
                cos_phi * cos_theta,
                sin_phi,
                cos_phi * sin_theta
            ])

            # Apply the combined rotation to the point and its normal
            P_rot = R @ P_orig
            N_rot = R @ N_orig
            
            # Normalize the rotated normal for lighting calculation
            N_rot_normalized = N_rot / np.linalg.norm(N_rot)

            # Calculate luminance
            luminance = np.dot(N_rot_normalized, light_source)
            
            # If the surface is visible (facing the light source)
            if luminance > 0:
                # Project the 3D point onto the 2D screen
                x_proj = int(screen_width / 2 + screen_scale * P_rot[0])
                y_proj = int(screen_height / 2 - screen_scale * P_rot[1] * 0.5) # Adjust for char aspect ratio

                if 0 <= y_proj < screen_height and 0 <= x_proj < screen_width:
                    # Use Z-buffer to handle occlusion
                    if P_rot[2] < z_buffer[y_proj, x_proj]:
                        z_buffer[y_proj, x_proj] = P_rot[2]
                        # Map luminance to a shade character
                        # High luminance -> low index -> '░'
                        # Low luminance -> high index -> '█'
                        shade_index = int((1 - luminance) * (len(shades) - 0.001))
                        output[y_proj, x_proj] = shades[shade_index]

    # --- Step 4: Print the final output ---
    print(f"Result for rotation (X={rx}, Y={ry}, Z={rz}):")
    print("-" * screen_width)
    for row in output:
        print("".join(row))
    print("-" * screen_width)

solve_torus_rotation()