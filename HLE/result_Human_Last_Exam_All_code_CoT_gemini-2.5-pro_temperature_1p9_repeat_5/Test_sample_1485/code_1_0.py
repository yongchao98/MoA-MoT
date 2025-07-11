import numpy as np

def render_rotated_torus():
    """
    This function models a torus, rotates it in 3D space, and renders
    an ASCII representation of the result from a fixed viewpoint.
    """
    # Torus parameters
    R = 2.0  # Major radius (center of torus to center of tube)
    r = 1.0  # Minor radius (radius of the tube)

    # Rotation angles in degrees. The problem states positive is clockwise.
    angle_x_deg = 140
    angle_y_deg = 75
    angle_z_deg = 35

    # Standard rotation matrices are counter-clockwise, so we negate the angles.
    # Convert degrees to radians for numpy's trigonometric functions.
    angle_x = np.radians(-angle_x_deg)
    angle_y = np.radians(-angle_y_deg)
    angle_z = np.radians(-angle_z_deg)

    # Define screen dimensions for the output ASCII art.
    width, height = 52, 22

    # A scaling factor to make the torus fit well on the screen.
    scale = 8

    # A correction factor for character aspect ratio (characters are taller than wide).
    aspect_ratio_correction = 0.5

    # Define the rotation matrix for the X-axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(angle_x), -np.sin(angle_x)],
        [0, np.sin(angle_x), np.cos(angle_x)]
    ])

    # Define the rotation matrix for the Y-axis
    Ry = np.array([
        [np.cos(angle_y), 0, np.sin(angle_y)],
        [0, 1, 0],
        [-np.sin(angle_y), 0, np.cos(angle_y)]
    ])

    # Define the rotation matrix for the Z-axis
    Rz = np.array([
        [np.cos(angle_z), -np.sin(angle_z), 0],
        [np.sin(angle_z), np.cos(angle_z), 0],
        [0, 0, 1]
    ])

    # Combine the matrices. Rotations are applied in order: X, then Y, then Z.
    R_matrix = Rz @ Ry @ Rx

    # Initialize the output grid (canvas) and a depth buffer (z-buffer).
    # The canvas stores the characters, and the z-buffer stores the depth of each point.
    canvas = np.full((height, width), ' ')
    z_buffer = np.full((height, width), np.inf)

    # Generate points on the torus surface by iterating over two angles, u and v.
    num_steps_u = 120
    num_steps_v = 120

    for u in np.linspace(0, 2 * np.pi, num_steps_u):
        for v in np.linspace(0, 2 * np.pi, num_steps_v):
            # Parametric equations for a torus with its major circle in the XZ-plane
            # This orients the hole along the Y-axis, matching the initial view.
            p = np.array([
                (R + r * np.cos(v)) * np.cos(u),
                r * np.sin(v),
                -(R + r * np.cos(v)) * np.sin(u)
            ])

            # Apply the combined rotation to the point
            p_rotated = R_matrix @ p
            x_rot, y_rot, z_rot = p_rotated[0], p_rotated[1], p_rotated[2]

            # Project the 3D point onto the 2D screen (orthographic projection)
            screen_x = int(width / 2 + scale * x_rot)
            screen_y = int(height / 2 - scale * y_rot * aspect_ratio_correction)

            # Z-buffering: If this point is on screen and closer than what's
            # already there, record its depth.
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if z_rot < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = z_rot

    # Map the depth values to the shading characters.
    shades = ['░', '▒', '▓', '█']  # Near (light) to far (dark)
    
    # Get the min and max depth values from the buffer to normalize them.
    finite_z = z_buffer[z_buffer != np.inf]
    if len(finite_z) > 0:
        min_z, max_z = finite_z.min(), finite_z.max()
        z_range = max_z - min_z
        if z_range == 0: z_range = 1 # Avoid division by zero if view is flat.

        # For each point on the canvas, calculate its shade based on its depth.
        for y in range(height):
            for x in range(width):
                if z_buffer[y, x] != np.inf:
                    depth = z_buffer[y, x]
                    # Normalize depth from 0 to 1
                    norm_depth = (depth - min_z) / z_range
                    # Choose shade index. The -0.001 is a trick to prevent index-out-of-bounds.
                    shade_index = int(norm_depth * (len(shades) - 0.001))
                    canvas[y, x] = shades[shade_index]

    # Print the equation and the final rendered view.
    print(f"Applying rotations X={angle_x_deg}, Y={angle_y_deg}, Z={angle_z_deg}:")
    print("\n".join("".join(row) for row in canvas))

if __name__ == '__main__':
    render_rotated_torus()