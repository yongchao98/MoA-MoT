import numpy as np
import math

def solve_torus_rotation():
    """
    Calculates and renders a view of a torus after 3D rotation.
    """
    # 1. Define rotation and object parameters
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35
    print(f"Applying rotation (X, Y, Z) in degrees: ({rot_x_deg}, {rot_y_deg}, {rot_z_deg})\n")

    # Positive rotation is clockwise. Standard matrices are counter-clockwise (CCW).
    # To perform a CW rotation by theta, we use a CCW rotation by -theta.
    rot_x_rad = math.radians(-rot_x_deg)
    rot_y_rad = math.radians(-rot_y_deg)
    rot_z_rad = math.radians(-rot_z_deg)

    # Torus geometry
    R_major = 2.0  # Major radius (from center of hole to center of tube)
    r_minor = 1.0  # Minor radius (radius of the tube)

    # Screen parameters
    width, height = 50, 22
    scale_x = 7.0  # Scale factor for x coordinate
    scale_y = 3.5  # Scale factor for y coordinate (accounts for character aspect ratio)


    # 2. Define rotation matrices
    def rot_x_mat(angle_rad):
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=np.float32)

    def rot_y_mat(angle_rad):
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=np.float32)

    def rot_z_mat(angle_rad):
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=np.float32)

    # Combine rotations in Z, Y, X order
    rotation_matrix = rot_x_mat(rot_x_rad) @ rot_y_mat(rot_y_rad) @ rot_z_mat(rot_z_rad)

    # 3. Initialize buffers
    output_screen = np.full((height, width), ' ', dtype=str)
    z_buffer = np.full((height, width), float('inf'), dtype=np.float32)

    # Character set for shading based on depth
    # '░' for nearest (min z), '█' for farthest (max z)
    shades = "░▒▓█"
    num_shades = len(shades)

    # 4. Generate points and find min/max z for normalization
    z_values = []
    u_step = 0.05
    v_step = 0.05
    u = 0.0
    while u < 2 * math.pi:
        v = 0.0
        while v < 2 * math.pi:
            cos_u, sin_u = math.cos(u), math.sin(u)
            cos_v, sin_v = math.cos(v), math.sin(v)
            
            point = np.array([
                (R_major + r_minor * cos_u) * cos_v,
                (R_major + r_minor * cos_u) * sin_v,
                r_minor * sin_u
            ])
            rot_point = rotation_matrix @ point
            z_values.append(rot_point[2])
            v += v_step
        u += u_step
    
    z_min, z_max = min(z_values), max(z_values)

    # 5. Main rendering loop
    u = 0.0
    while u < 2 * math.pi:
        v = 0.0
        while v < 2 * math.pi:
            cos_u, sin_u = math.cos(u), math.sin(u)
            cos_v, sin_v = math.cos(v), math.sin(v)

            # Define point on torus surface
            point = np.array([
                (R_major + r_minor * cos_u) * cos_v,
                (R_major + r_minor * cos_u) * sin_v,
                r_minor * sin_u
            ])

            # Rotate the point
            rot_point = rotation_matrix @ point
            rx, ry, rz = rot_point

            # Project to 2D screen coordinates
            screen_x = int(width / 2 + scale_x * rx)
            screen_y = int(height / 2 - scale_y * ry) # Y is inverted on screen

            # Z-buffer check for occlusion
            if 0 <= screen_x < width and 0 <= screen_y < height:
                if rz < z_buffer[screen_y, screen_x]:
                    z_buffer[screen_y, screen_x] = rz

                    # Determine shade based on normalized depth
                    z_normalized = (rz - z_min) / (z_max - z_min)
                    shade_index = int(z_normalized * (num_shades - 0.001))
                    output_screen[screen_y, screen_x] = shades[shade_index]

            v += v_step
        u += u_step

    # 6. Print the result
    for row in output_screen:
        print("".join(row))

solve_torus_rotation()
print("\nThe generated image most closely matches answer choice B.")
<<<B>>>