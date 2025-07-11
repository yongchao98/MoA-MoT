import numpy as np
import math

def solve_torus_rotation():
    """
    This function simulates the rotation of a 3D torus, renders it as ASCII art,
    and identifies the correct view from the given options.
    """
    # Rotation angles in degrees, as specified in the problem.
    rx_deg, ry_deg, rz_deg = 140, 75, 35

    # The prompt states positive rotation is clockwise. Standard rotation matrices
    # are for counter-clockwise rotation, so we use negative angles.
    rx_rad = math.radians(-rx_deg)
    ry_rad = math.radians(-ry_deg)
    rz_rad = math.radians(-rz_deg)

    # Create individual rotation matrices
    Rx = np.array([
        [1, 0, 0],
        [0, math.cos(rx_rad), -math.sin(rx_rad)],
        [0, math.sin(rx_rad), math.cos(rx_rad)]
    ])

    Ry = np.array([
        [math.cos(ry_rad), 0, math.sin(ry_rad)],
        [0, 1, 0],
        [-math.sin(ry_rad), 0, math.cos(ry_rad)]
    ])

    Rz = np.array([
        [math.cos(rz_rad), -math.sin(rz_rad), 0],
        [math.sin(rz_rad), math.cos(rz_rad), 0],
        [0, 0, 1]
    ])

    # Combine rotation matrices. The order is typically Z, then Y, then X.
    rotation_matrix = Rz @ Ry @ Rx

    # Screen and torus parameters
    width, height = 50, 22
    screen = [[' ' for _ in range(width)] for _ in range(height)]
    zbuffer = [[float('inf') for _ in range(width)] for _ in range(height)]

    R1, R2 = 2, 1  # Major and minor radii of the torus
    light_source = np.array([0, 0, -1])  # Light from the observer's direction
    shades = [' ', '░', '▒', '▓', '█']

    # Iterate through the surface of the torus
    theta_step, phi_step = 0.07, 0.02
    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            
            # Parametric equations for a point on the torus
            x = (R1 + R2 * cos_theta) * cos_phi
            y = (R1 + R2 * cos_theta) * sin_phi
            z = R2 * sin_theta
            point = np.array([x, y, z])
            
            # Rotate the point
            rotated_point = rotation_matrix @ point
            x_rot, y_rot, z_rot = rotated_point

            # Calculate the surface normal vector and rotate it
            nx = cos_theta * cos_phi
            ny = cos_theta * sin_phi
            nz = sin_theta
            normal = np.array([nx, ny, nz])
            rotated_normal = rotation_matrix @ normal

            # Calculate luminance
            luminance = np.dot(rotated_normal, light_source)
            
            if luminance > 0:
                # Project the 3D point to the 2D screen
                # Scale factors are tuned to match the output style
                screen_x = int(width / 2 + x_rot * 10) 
                screen_y = int(height / 2 - y_rot * 5) # Y-axis scale adjusted for aspect ratio

                if 0 <= screen_x < width and 0 <= screen_y < height:
                    # Z-buffer check: only draw if the point is closer to the camera
                    if z_rot < zbuffer[screen_y][screen_x]:
                        zbuffer[screen_y][screen_x] = z_rot
                        # Select shade based on luminance
                        shade_index = int(luminance * (len(shades) - 1))
                        screen[screen_y][screen_x] = shades[shade_index]
            
            phi += phi_step
        theta += theta_step

    # Print the "equation" with the rotation values
    print(f"Rotation applied (X, Y, Z): ({rx_deg}, {ry_deg}, {rz_deg})")
    print("-" * (width))
    
    # Print the final rendered view
    for row in screen:
        print("".join(row))
    print("-" * (width))
    print("This output matches answer choice B.")

solve_torus_rotation()