import numpy as np
import math

def solve():
    """
    This function simulates the rotation of a 3D torus and renders it as ASCII art.
    """
    # 1. Define parameters for the torus, screen, and rotations.
    R = 7.0  # Major radius of the torus
    r = 3.0  # Minor radius of the torus

    # Screen dimensions in characters
    width, height = 45, 22

    # Rotation angles from the problem description
    angle_x_deg = 140
    angle_y_deg = 75
    angle_z_deg = 35

    print(f"Applying rotation (X, Y, Z) in degrees: ({angle_x_deg}, {angle_y_deg}, {angle_z_deg})")
    print("-" * 30)

    # Convert angles to radians for trigonometric functions
    angle_x = math.radians(angle_x_deg)
    angle_y = math.radians(angle_y_deg)
    angle_z = math.radians(angle_z_deg)

    # 2. Define helper functions for clockwise rotation matrices.
    def get_rotation_x_cw(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])

    def get_rotation_y_cw(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])

    def get_rotation_z_cw(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

    # Combine rotation matrices. The order is X, then Y, then Z.
    # A point P is transformed as P' = Rz * Ry * Rx * P
    rotation_matrix = get_rotation_z_cw(angle_z) @ get_rotation_y_cw(angle_y) @ get_rotation_x_cw(angle_x)

    # 3. Initialize buffers for projection.
    # z_buffer stores the depth of the closest point for each screen character.
    # char_buffer stores the character to be displayed.
    z_buffer = np.full((height, width), float('inf'))
    char_buffer = np.full((height, width), ' ', dtype=str)

    # Define point density for a smooth surface
    num_steps_u = 150  # Steps around the major circle
    num_steps_v = 80   # Steps around the minor circle (the tube)

    # 4. Generate, rotate, and project points onto the screen.
    for u in np.linspace(0, 2 * math.pi, num_steps_u):
        for v in np.linspace(0, 2 * math.pi, num_steps_v):
            # Parametric equations for a torus with its hole along the Y-axis
            x = (R + r * math.cos(v)) * math.cos(u)
            y = r * math.sin(v)
            z = (R + r * math.cos(v)) * math.sin(u)
            point = np.array([x, y, z])

            # Apply the combined rotation
            rotated_point = rotation_matrix @ point
            rx, ry, rz = rotated_point[0], rotated_point[1], rotated_point[2]

            # Project onto a 2D screen (orthographic projection)
            # Scale and shift to fit the character grid
            scale = 1.9
            screen_x = int(width / 2 + rx * scale)
            screen_y = int(height / 2 - ry * scale)  # Y-axis is inverted in screen coordinates

            # Z-buffering: check if the point is visible and closer than existing points
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if rz < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = rz

    # 5. Map depth values to shading characters.
    shades = ['░', '▒', '▓', '█']
    # Find the min and max z values of visible points to normalize the depth
    min_z = np.min(z_buffer[z_buffer != float('inf')])
    max_z = np.max(z_buffer[z_buffer != float('inf')])
    z_range = max_z - min_z if max_z > min_z else 1.0

    for y in range(height):
        for x in range(width):
            if z_buffer[y, x] != float('inf'):
                # Normalize depth to a value between 0 and 1
                normalized_z = (z_buffer[y, x] - min_z) / z_range
                # Select the character based on the normalized depth
                shade_index = int(normalized_z * (len(shades) - 0.001))
                char_buffer[y, x] = shades[shade_index]

    # 6. Print the final ASCII art.
    for row in char_buffer:
        print("".join(row).rstrip())

solve()