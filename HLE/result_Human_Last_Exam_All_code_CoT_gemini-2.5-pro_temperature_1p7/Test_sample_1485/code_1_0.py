import math
import numpy as np

def solve_torus_rotation():
    """
    This function generates an ASCII view of a torus after specified rotations
    and prints the result to the console.
    """
    
    # Rotation angles in degrees from the problem description
    angle_x, angle_y, angle_z = 140, 75, 35

    print(f"Calculating view for rotation (X, Y, Z) in degrees: ({angle_x}, {angle_y}, {angle_z})")
    print("-" * 60)

    # Helper function for creating clockwise rotation matrices
    def get_rotation_matrix(angle_deg, axis):
        theta = math.radians(angle_deg)
        c, s = math.cos(theta), math.sin(theta)
        if axis == 'x':
            return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
        elif axis == 'y':
            return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])
        elif axis == 'z':
            return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        return np.identity(3)

    # Torus and screen parameters
    R, r = 2.0, 1.0  # Major and minor radii
    width, height = 52, 22
    screen_buffer = [[' ' for _ in range(width)] for _ in range(height)]
    z_buffer = [[-float('inf') for _ in range(width)] for _ in range(height)]

    # Combine rotation matrices in the specified order: X, then Y, then Z
    Rx = get_rotation_matrix(angle_x, 'x')
    Ry = get_rotation_matrix(angle_y, 'y')
    Rz = get_rotation_matrix(angle_z, 'z')
    M = Rz @ Ry @ Rx

    # Light/View vector (from observer at +Z)
    light_vector = np.array([0, 0, 1])

    # Shading characters based on luminance (from low to high)
    shades = "░▒▓█"

    # Iterate through the torus surface's parameters (u and v)
    u_steps, v_steps = 100, 100
    for i in range(u_steps):
        u = 2 * math.pi * i / u_steps
        cos_u, sin_u = math.cos(u), math.sin(u)
        for j in range(v_steps):
            v = 2 * math.pi * j / v_steps
            cos_v, sin_v = math.cos(v), math.sin(v)

            # Define point P and normal N for a torus lying on the XY plane
            P = np.array([
                (R + r * cos_u) * cos_v,
                (R + r * cos_u) * sin_v,
                r * sin_u
            ])
            N = np.array([
                cos_u * cos_v,
                cos_u * sin_v,
                sin_u
            ])

            # Apply the rotation to the point and its normal vector
            P_rot = M @ P
            N_rot = M @ N

            # Calculate luminance and check if the surface faces the camera
            luminance = np.dot(N_rot, light_vector)
            if luminance > 0:
                x_proj, y_proj, z_proj = P_rot[0], P_rot[1], P_rot[2]

                # Project to 2D screen coordinates and scale
                # The scale is tweaked to best match the answer choices' size.
                scale = height / (2 * (R + r)) * 1.6
                # Adjust for character aspect ratio (typically 2:1 height:width)
                sx = int(width / 2 + scale * 2.0 * x_proj)
                sy = int(height / 2 - scale * y_proj)
                
                # Update buffers if the point is visible and closer than previous points
                if 0 <= sx < width and 0 <= sy < height:
                    if z_proj > z_buffer[sy][sx]:
                        z_buffer[sy][sx] = z_proj
                        shade_index = min(int(luminance * len(shades)), len(shades) - 1)
                        screen_buffer[sy][sx] = shades[shade_index]

    # Print the final rendered view
    # The output includes blank lines at top and bottom to match the answer's format
    print()
    for row in screen_buffer:
        print("".join(row))
    print()


solve_torus_rotation()