import numpy as np

def solve_torus_rotation():
    """
    This function generates and rotates a 3D torus, then renders it to the console
    using ASCII characters based on depth, matching the specified rotation.
    """
    # Screen and torus parameters
    width, height = 50, 22
    R, r = 2.0, 1.0  # Major and minor radii of the torus
    scale = 8.5
    # Adjusted offset to match the position of the answer choice
    offset_x = width / 2 + 2 
    offset_y = height / 2

    # --- The Rotation Equation ---
    # The final rotation is a composition of rotations around X, Y, and Z axes.
    # P_rotated = Rz(gamma) * Ry(beta) * Rx(alpha) * P_original
    # Where P is a point on the torus, and Rx, Ry, Rz are the rotation matrices.
    x_rot_deg, y_rot_deg, z_rot_deg = 140, 75, 35
    
    print("Generating view for a torus after rotation.")
    print(f"Rotation around X-axis (alpha): {x_rot_deg} degrees")
    print(f"Rotation around Y-axis (beta): {y_rot_deg} degrees")
    print(f"Rotation around Z-axis (gamma): {z_rot_deg} degrees")
    print("-" * (width+2))

    # Convert angles to radians. Negate for clockwise rotation using standard
    # counter-clockwise rotation matrices.
    ax = -np.radians(x_rot_deg)
    ay = -np.radians(y_rot_deg)
    az = -np.radians(z_rot_deg)

    # Create individual rotation matrices
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(ax), -np.sin(ax)],
        [0, np.sin(ax), np.cos(ax)]
    ])
    Ry = np.array([
        [np.cos(ay), 0, np.sin(ay)],
        [0, 1, 0],
        [-np.sin(ay), 0, np.cos(ay)]
    ])
    Rz = np.array([
        [np.cos(az), -np.sin(az), 0],
        [np.sin(az), np.cos(az), 0],
        [0, 0, 1]
    ])
    
    # Combine matrices (rotations applied in order X, then Y, then Z)
    rotation_matrix = Rz @ Ry @ Rx

    # Initialize buffers for rendering
    z_buffer = np.full((height, width), np.inf)
    output_buffer = np.full((height, width), ' ')

    # Generate points on the torus surface
    theta_spacing = 0.07
    phi_spacing = 0.02
    for theta in np.arange(0, 2 * np.pi, theta_spacing):
        for phi in np.arange(0, 2 * np.pi, phi_spacing):
            # Original point on the torus
            p_original = np.array([
                (R + r * np.cos(theta)) * np.cos(phi),
                (R + r * np.cos(theta)) * np.sin(phi),
                r * np.sin(theta)
            ])
            
            # Apply rotation
            p_rotated = rotation_matrix @ p_original
            x, y, z = p_rotated[0], p_rotated[1], p_rotated[2]
            
            # Project to 2D screen coordinates
            # Y is inverted because screen coordinates have y=0 at the top
            screen_x = int(offset_x + scale * x)
            screen_y = int(offset_y - scale * y)

            # Update z-buffer if the point is visible and closer
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if z < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = z

    # Determine depth range for shading
    finite_z = z_buffer[z_buffer != np.inf]
    if len(finite_z) == 0:
        return # Avoid errors if nothing is rendered
        
    min_z, max_z = finite_z.min(), finite_z.max()
    depth_range = max_z - min_z
    
    # Define shading characters from nearest (light) to farthest (dark)
    chars = ['░', '▒', '▓', '█']
    
    # Fill the output buffer based on depth
    for r_idx in range(height):
        for c_idx in range(width):
            if z_buffer[r_idx, c_idx] != np.inf:
                # Normalize depth to a [0, 1] range
                normalized_z = (z_buffer[r_idx, c_idx] - min_z) / depth_range if depth_range > 0 else 0
                # Map normalized depth to a character index
                char_index = int(normalized_z * (len(chars) - 0.001))
                output_buffer[r_idx, c_idx] = chars[char_index]

    # Print the final rendered image
    for row in output_buffer:
        print("".join(row))

solve_torus_rotation()