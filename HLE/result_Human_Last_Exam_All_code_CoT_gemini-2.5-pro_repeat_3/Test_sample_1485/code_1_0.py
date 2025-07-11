import math
import numpy as np

def solve_torus_rotation():
    """
    This function models, rotates, and renders a torus in ASCII art
    based on the specified rotation angles.
    """
    
    # --- Configuration ---
    # Screen dimensions for the ASCII art
    width, height = 44, 22
    
    # Torus geometry: R is the major radius, r is the minor radius
    R, r = 2.0, 1.0
    
    # Rotation angles in degrees as specified in the problem
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35
    
    # Convert angles to radians. The problem states positive rotation is
    # clockwise, while standard math functions use counter-clockwise.
    # Therefore, we use negative angles for the rotation.
    rot_x_rad = math.radians(-rot_x_deg)
    rot_y_rad = math.radians(-rot_y_deg)
    rot_z_rad = math.radians(-rot_z_deg)
    
    # Projection settings
    screen_scale = 8  # Scales the torus to fit the screen
    y_aspect = 0.5    # Corrects for the aspect ratio of text characters
    
    # --- Buffers ---
    # z_buffer stores the depth of the closest point for each pixel
    z_buffer = np.full((height, width), float('inf'))
    # output_buffer will store the final characters for rendering
    output_buffer = np.full((height, width), ' ')
    
    # --- Rotation and Projection Loop ---
    # Iterate through the torus surface using angles u and v
    for u in np.linspace(0, 2 * math.pi, 200):  # Major circle angle
        for v in np.linspace(0, 2 * math.pi, 200):  # Minor circle angle
            
            # --- 1. Model Point on Torus ---
            p_x = (R + r * math.cos(v)) * math.cos(u)
            p_y = (R + r * math.cos(v)) * math.sin(u)
            p_z = r * math.sin(v)
            
            # --- 2. Apply 3D Rotations ---
            # Rotate around X-axis
            y_rot1 = p_y * math.cos(rot_x_rad) - p_z * math.sin(rot_x_rad)
            z_rot1 = p_y * math.sin(rot_x_rad) + p_z * math.cos(rot_x_rad)
            x_rot1 = p_x
            
            # Rotate around Y-axis
            x_rot2 = x_rot1 * math.cos(rot_y_rad) + z_rot1 * math.sin(rot_y_rad)
            z_rot2 = -x_rot1 * math.sin(rot_y_rad) + z_rot1 * math.cos(rot_y_rad)
            y_rot2 = y_rot1
            
            # Rotate around Z-axis
            x_final = x_rot2 * math.cos(rot_z_rad) - y_rot2 * math.sin(rot_z_rad)
            y_final = x_rot2 * math.sin(rot_z_rad) + y_rot2 * math.cos(rot_z_rad)
            z_final = z_rot2
            
            # --- 3. Project to 2D Screen ---
            screen_x = int(width / 2 + screen_scale * x_final)
            screen_y = int(height / 2 - screen_scale * y_aspect * y_final)
            
            # --- 4. Handle Occlusion (Z-buffering) ---
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if z_final < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = z_final

    # --- 5. Render with Depth Shading ---
    # Find the min and max depth values of the visible points
    min_z = np.min(z_buffer[z_buffer != float('inf')])
    max_z = np.max(z_buffer)
    depth_range = max_z - min_z
    if depth_range == 0: depth_range = 1 # Avoid division by zero
        
    depth_shades = ['░', '▒', '▓', '█']

    for y in range(height):
        for x in range(width):
            if z_buffer[y, x] != float('inf'):
                # Normalize depth: 0 for nearest (min_z), 1 for farthest (max_z)
                normalized_depth = (z_buffer[y, x] - min_z) / depth_range
                # Select shade based on normalized depth
                shade_index = int(normalized_depth * (len(depth_shades) - 0.001))
                output_buffer[y, x] = depth_shades[shade_index]

    # --- Print Final Result ---
    print(f"Generating view for rotation (X, Y, Z) = ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees.\n")

    # Trim empty border lines for cleaner output
    non_empty_rows = [i for i, row in enumerate(output_buffer) if not np.all(row == ' ')]
    if not non_empty_rows:
        return
        
    for i in range(min(non_empty_rows), max(non_empty_rows) + 1):
        print("".join(output_buffer[i]))

solve_torus_rotation()