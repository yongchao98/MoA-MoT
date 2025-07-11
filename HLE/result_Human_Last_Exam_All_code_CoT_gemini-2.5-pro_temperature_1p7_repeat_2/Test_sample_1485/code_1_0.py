import numpy as np
import math

def solve():
    """
    This function simulates the rotation of a torus and renders the result.
    """
    # --- Parameters ---
    SCREEN_WIDTH = 48
    SCREEN_HEIGHT = 22
    R1 = 10  # Major radius of the torus
    R2 = 4   # Minor radius of the torus
    SCALE_X = 2.0  # Correction for character aspect ratio in terminals
    SCALE_Y = 1.0

    # --- Rotation Angles (in degrees) ---
    angle_x_deg = 140
    angle_y_deg = 75
    angle_z_deg = 35

    ax_rad, ay_rad, az_rad = map(math.radians, [angle_x_deg, angle_y_deg, angle_z_deg])

    # --- Rotation Matrices (Clockwise) ---
    cx, sx = math.cos(ax_rad), math.sin(ax_rad)
    Rx = np.array([[1, 0, 0], [0, cx, sx], [0, -sx, cx]])

    cy, sy = math.cos(ay_rad), math.sin(ay_rad)
    Ry = np.array([[cy, 0, -sy], [0, 1, 0], [sy, 0, cy]])

    cz, sz = math.cos(az_rad), math.sin(az_rad)
    Rz = np.array([[cz, sz, 0], [-sz, cz, 0], [0, 0, 1]])

    # Combine rotations (extrinsic order: Z then Y then X)
    R = Rz @ Ry @ Rx

    # --- Buffers ---
    # z_buffer stores the depth of the closest point at each screen cell
    z_buffer = np.full((SCREEN_HEIGHT, SCREEN_WIDTH), float('inf'))
    # output_buffer stores the final character for each screen cell
    output_buffer = np.full((SCREEN_HEIGHT, SCREEN_WIDTH), ' ')

    # --- Torus Point Generation and Rotation ---
    theta_step = 0.07  # Step size for torus main circle
    phi_step = 0.035   # Step size for torus tube circle

    for theta in np.arange(0, 2 * math.pi, theta_step):
        for phi in np.arange(0, 2 * math.pi, phi_step):
            # 1. Point on a standard torus (hole along z-axis)
            x_std = (R1 + R2 * math.cos(phi)) * math.cos(theta)
            y_std = (R1 + R2 * math.cos(phi)) * math.sin(theta)
            z_std = R2 * math.sin(phi)
            
            # 2. Rotate to initial orientation (hole along y-axis)
            # This is a +90 degree CCW rotation around the x-axis
            x_init = x_std
            y_init = -z_std
            z_init = y_std
            p_initial = np.array([x_init, y_init, z_init])

            # 3. Apply the main user-specified rotation
            p_rotated = R @ p_initial
            x, y, z = p_rotated[0], p_rotated[1], p_rotated[2]

            # 4. Project onto 2D screen (Orthographic Projection)
            screen_scale = min(SCREEN_WIDTH, SCREEN_HEIGHT) / (R1 + R2) / 2.5
            x_proj = int(SCREEN_WIDTH / 2 + SCALE_X * screen_scale * x)
            # Y-axis is inverted for screen coordinates (0 is top)
            y_proj = int(SCREEN_HEIGHT / 2 - SCALE_Y * screen_scale * y)

            # 5. Z-buffer check
            if 0 <= y_proj < SCREEN_HEIGHT and 0 <= x_proj < SCREEN_WIDTH:
                if z < z_buffer[y_proj, x_proj]:
                    z_buffer[y_proj, x_proj] = z

    # --- Shading Pass ---
    # Find the range of z-depths for visible points
    visible_zs = z_buffer[z_buffer != float('inf')]
    if len(visible_zs) > 0:
        min_z, max_z = np.min(visible_zs), np.max(visible_zs)
        z_range = max_z - min_z if max_z > min_z else 1.0

        shades = "░▒▓█"
        num_shades = len(shades)
        
        for r in range(SCREEN_HEIGHT):
            for c in range(SCREEN_WIDTH):
                if z_buffer[r, c] != float('inf'):
                    # Lighter shades for nearer points (smaller z)
                    depth = z_buffer[r, c]
                    shade_index = int(((depth - min_z) / z_range) * (num_shades - 0.001))
                    output_buffer[r, c] = shades[shade_index]

    # --- Print Final Output ---
    for row in output_buffer:
        print("".join(list(row)))
        
    # As requested, output each number in the final equation
    print("\nRotation equation values used:")
    print(f"X-axis rotation: {angle_x_deg} degrees")
    print(f"Y-axis rotation: {angle_y_deg} degrees")
    print(f"Z-axis rotation: {angle_z_deg} degrees")


solve()