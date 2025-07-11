import numpy as np
import math

def generate_rotated_torus_view():
    """
    This function simulates the rotation of a 3D torus and renders a 2D ASCII view.
    """
    # Screen dimensions for the output ASCII art
    screen_width = 50
    screen_height = 21

    # Torus parameters: R1 is the major radius, R2 is the minor radius.
    R1 = 2.0
    R2 = 1.0

    # Observer's distance from the origin. Affects perspective.
    K2 = 5.0
    # A scaling factor for the projection.
    K1 = screen_width * K2 * 3 / (8 * (R1 + R2))

    # Rotation angles in degrees as specified in the problem
    angle_x_deg, angle_y_deg, angle_z_deg = 140, 75, 35

    # Convert degrees to radians for trigonometric functions
    angle_x = math.radians(angle_x_deg)
    angle_y = math.radians(angle_y_deg)
    angle_z = math.radians(angle_z_deg)

    # To perform clockwise rotations using standard counter-clockwise rotation matrices,
    # we use the negative of the angles.
    sin_x, cos_x = math.sin(-angle_x), math.cos(-angle_x)
    sin_y, cos_y = math.sin(-angle_y), math.cos(-angle_y)
    sin_z, cos_z = math.sin(-angle_z), math.cos(-angle_z)

    # Define the rotation matrices for each axis
    Rx = np.array([[1, 0, 0], [0, cos_x, -sin_x], [0, sin_x, cos_x]])
    Ry = np.array([[cos_y, 0, sin_y], [0, 1, 0], [-sin_y, 0, cos_y]])
    Rz = np.array([[cos_z, -sin_z, 0], [sin_z, cos_z, 0], [0, 0, 1]])

    # Combine the rotation matrices. Rotations are applied in the order: X, then Y, then Z.
    R_combined = Rz @ Ry @ Rx
    
    # Define a light source. To match the visual style of "nearer is lighter",
    # we use a light source pointing from the viewer towards the object.
    # In a left-handed system (Z-in), this is (0, 0, -1).
    light_source = np.array([0, 0, -1])
    light_source = light_source / np.linalg.norm(light_source)

    # Initialize buffers for rendering
    z_buffer = np.full((screen_height, screen_width), -np.inf)
    char_buffer = np.full((screen_height, screen_width), ' ', dtype='<U1')

    # Iterate through the torus's surface using parametric angles theta and phi
    for theta in np.arange(0, 2 * math.pi, 0.04):
        for phi in np.arange(0, 2 * math.pi, 0.02):
            # Parametric equations for the initial torus (hole along Y-axis)
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)

            # Initial point coordinates
            x0 = (R1 + R2 * cos_theta) * cos_phi
            y0 = R2 * sin_theta
            z0 = (R1 + R2 * cos_theta) * sin_phi
            p0 = np.array([x0, y0, z0])

            # Initial surface normal vector
            nx0 = cos_theta * cos_phi
            ny0 = sin_theta
            nz0 = cos_theta * sin_phi
            n0 = np.array([nx0, ny0, nz0])
            
            # Apply the combined rotation to the point and the normal vector
            p_rotated = R_combined @ p0
            n_rotated = R_combined @ n0
            x, y, z = p_rotated

            # Perspective projection
            z_proj = z + K2
            if z_proj == 0: continue
            
            ooz = 1.0 / z_proj # 1/z is used for the depth buffer

            # Project to 2D screen coordinates
            # A 0.5 factor on y-coord helps adjust for non-square character aspect ratio
            xp = int(screen_width / 2 + K1 * x * ooz)
            yp = int(screen_height / 2 - K1 * y * ooz * 0.5)

            # Z-buffering: check if the point is on screen and closer than what's already there
            if 0 <= xp < screen_width and 0 <= yp < screen_height:
                if ooz > z_buffer[yp, xp]:
                    z_buffer[yp, xp] = ooz
                    
                    # Calculate luminance from the dot product of the light source and the rotated normal
                    luminance = np.dot(n_rotated, light_source)
                    
                    # Choose character based on luminance
                    shading_chars = '░▒▓█'
                    if luminance > 0:
                        char_index = int(luminance * (len(shading_chars) - 0.001))
                        char_buffer[yp, xp] = shading_chars[char_index]
                    else:
                        # Surfaces facing away from the light get the 'nearest' char to simulate outline
                        char_buffer[yp, xp] = shading_chars[0]

    # Print the resulting ASCII art
    print("--- Simulated View ---")
    for row in char_buffer:
        print("".join(row))
    
    # As requested, print the numbers used in the calculation
    print("\n--- Rotation Parameters ---")
    print(f"Rotation around X-axis: {angle_x_deg} degrees")
    print(f"Rotation around Y-axis: {angle_y_deg} degrees")
    print(f"Rotation around Z-axis: {angle_z_deg} degrees")

# Run the simulation and display the result
generate_rotated_torus_view()