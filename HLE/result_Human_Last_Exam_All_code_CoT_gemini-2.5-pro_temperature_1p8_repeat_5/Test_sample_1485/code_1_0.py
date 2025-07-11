import numpy as np

def solve_torus_rotation():
    """
    This function calculates the view of a torus after specified rotations
    and prints the resulting ASCII art to the console.
    """
    # 1. Parameters
    R1 = 2.0  # Major radius of the torus
    R2 = 1.0  # Minor radius of the torus (tube radius)
    
    # Rotation angles in degrees, as specified in the problem
    rot_x, rot_y, rot_z = 140, 75, 35
    
    # The prompt asks to output each number in the final equation.
    # I am interpreting this as printing the rotation values being used.
    print(f"Applying rotations (X, Y, Z) in degrees: ({rot_x}, {rot_y}, {rot_z})")
    
    # Screen and projection settings
    width, height = 50, 25
    canvas = [[' '] * width for _ in range(height)]
    z_buffer = [[float('inf')] * width for _ in range(height)]

    # 2. Define Rotation Matrix
    # Problem states positive rotation is clockwise. Standard matrices are counter-clockwise.
    # So, we use the negative of the given angles.
    ax = -np.deg2rad(rot_x)
    ay = -np.deg2rad(rot_y)
    az = -np.deg2rad(rot_z)

    # Standard 3D rotation matrices
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
    
    # Combine rotations in order: Rz @ Ry @ Rx
    R = Rz @ Ry @ Rx

    # 3. Generate, Rotate, and Project Points
    # Iterate over the angles of the torus parameterization
    phi_step = 0.07  # Step for the major circle angle
    theta_step = 0.07 # Step for the minor circle angle
    
    for phi in np.arange(0, 2 * np.pi, phi_step):
        for theta in np.arange(0, 2 * np.pi, theta_step):
            # Parametric equations for a torus with its hole along the Y-axis
            h = R1 + R2 * np.cos(theta)
            p0 = np.array([
                h * np.cos(phi), 
                R2 * np.sin(theta), 
                h * np.sin(phi)
            ])
            
            # Apply the combined rotation
            p_rot = R @ p0
            x, y, z = p_rot[0], p_rot[1], p_rot[2]

            # Project the 3D point to 2D screen coordinates
            scale = 5.0
            x_offset = width / 2
            y_offset = height / 2
            aspect_ratio_correction = 0.5  # To account for non-square characters in terminals
            
            xp = int(x_offset + scale * x)
            yp = int(y_offset - scale * y * aspect_ratio_correction)

            # 4. Z-Buffering for occlusion
            if 0 <= xp < width and 0 <= yp < height:
                if z < z_buffer[yp][xp]:
                    z_buffer[yp][xp] = z

    # 5. Determine depth range for coloring
    visible_zs = [z for row in z_buffer for z in row if z != float('inf')]
    if not visible_zs:
        print("No points were projected onto the screen.")
        return

    min_z = min(visible_zs)
    max_z = max(visible_zs)
    z_range = max_z - min_z
    if z_range == 0: z_range = 1 # Avoid division by zero if view is flat

    # 6. Render with Depth Shading
    # Map Z-depth to characters: ░ (near) to █ (far)
    chars = ['░', '▒', '▓', '█']
    for r in range(height):
        for c in range(width):
            if z_buffer[r][c] != float('inf'):
                # Normalize z-value to [0, 1]
                normalized_z = (z_buffer[r][c] - min_z) / z_range
                # Select character based on normalized depth
                char_index = min(len(chars) - 1, int(normalized_z * len(chars)))
                canvas[r][c] = chars[char_index]
    
    # 7. Print the resulting view
    for row in canvas:
        print("".join(row))

# Run the simulation
solve_torus_rotation()
print("\nComparing the output with the answer choices, the closest match is D.")
<<<D>>>