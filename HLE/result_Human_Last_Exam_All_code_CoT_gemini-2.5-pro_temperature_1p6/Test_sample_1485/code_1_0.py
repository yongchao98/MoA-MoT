import numpy as np

def solve_torus_rotation():
    """
    This function generates and displays a rotated torus using ASCII art.
    It models a torus in 3D, applies specified rotations, projects it to a 2D
    plane, and uses a z-buffer for occlusion and depth-based shading.
    """
    # --- Parameters ---
    # Screen dimensions
    width, height = 50, 24
    
    # Torus geometry
    R, r = 2.2, 1.1  # Major and minor radii
    
    # Input rotation values in degrees as per the problem
    rot_x_deg_in, rot_y_deg_in, rot_z_deg_in = 140, 75, 35
    
    # Print the equation/input values
    print(f"Generating view for rotation (X, Y, Z): ({rot_x_deg_in}, {rot_y_deg_in}, {rot_z_deg_in})")
    
    # Note: "positive rotation turns clockwise" is specified.
    # Standard rotation matrices are counter-clockwise.
    # To get clockwise, we use the negative of the angle.
    rot_x_deg = -rot_x_deg_in
    rot_y_deg = -rot_y_deg_in
    rot_z_deg = -rot_z_deg_in

    # --- Math Setup ---
    # Screen and Z-buffer initialization
    z_buffer = np.full((height, width), float('inf'))
    output_buffer = np.full((height, width), ' ')
    
    # Convert angles to radians
    rot_x = np.radians(rot_x_deg)
    rot_y = np.radians(rot_y_deg)
    rot_z = np.radians(rot_z_deg)
    
    # Create rotation matrices
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
    
    # Combined rotation matrix for efficiency (Z*Y*X)
    R_total = Rz @ Ry @ Rx

    # --- Torus Generation and Projection ---
    theta_step = 0.05  # Steps for iterating around the tube
    phi_step = 0.05    # Steps for iterating around the main ring

    for theta in np.arange(0, 2 * np.pi, theta_step):
        for phi in np.arange(0, 2 * np.pi, phi_step):
            # 3D point on the original torus surface
            x = (R + r * np.cos(theta)) * np.cos(phi)
            y = (R + r * np.cos(theta)) * np.sin(phi)
            z = r * np.sin(theta)
            point = np.array([x, y, z])
            
            # Apply the combined rotation
            rotated_point = R_total @ point
            rx, ry, rz = rotated_point
            
            # Orthographic projection and scaling to screen coordinates
            scale = 5.5
            # Y-axis is bottom-to-top, so we flip y for screen display
            screen_x = int(width / 2 + scale * rx)
            screen_y = int(height / 2 - scale * ry)
            
            # Z-buffer check for occlusion
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if rz < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = rz

    # --- Shading Pass ---
    # Shading characters: Near (light) -> Far (dark)
    shades = '░▒▓█'
    
    min_z = np.min(z_buffer[z_buffer != float('inf')])
    max_z = np.max(z_buffer[z_buffer != float('inf')])
    
    if max_z == min_z: # Avoid division by zero if view is flat
        max_z += 1

    for y in range(height):
        for x in range(width):
            if z_buffer[y, x] != float('inf'):
                # Normalize z-depth from 0 (near) to 1 (far)
                norm_z = (z_buffer[y, x] - min_z) / (max_z - min_z)
                # Select character based on depth
                char_index = int(norm_z * (len(shades) - 0.001))
                output_buffer[y, x] = shades[char_index]

    # --- Print Final ASCII Art ---
    print("\n--- Generated View ---")
    for row in output_buffer:
        print("".join(row))
    print("----------------------")

# Execute the function to get the result
solve_torus_rotation()