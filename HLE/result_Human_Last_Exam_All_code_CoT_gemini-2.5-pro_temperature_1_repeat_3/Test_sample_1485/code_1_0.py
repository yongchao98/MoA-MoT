import numpy as np
import math

def solve_torus_rotation():
    """
    This function generates and displays an ASCII art representation of a torus
    after applying a specific 3D rotation.

    The plan is as follows:
    1.  Define the geometric and display parameters: major radius (R1), minor radius (R2) of the torus,
        and the dimensions of the output character grid (screen).
    2.  Define the rotation angles (140, 75, 35) for the X, Y, and Z axes respectively.
    3.  Convert the rotation angles from degrees to radians. A positive rotation is defined as clockwise,
        which is opposite to the standard mathematical convention. Therefore, we will use the negative of the angles.
    4.  Construct the individual rotation matrices for each axis (Rx, Ry, Rz) using the calculated radians.
    5.  Combine these matrices by multiplication (Rz @ Ry @ Rx) to get a single transformation matrix
        that applies the rotations in the specified order: X, then Y, then Z.
    6.  Initialize a 2D character buffer for the final output and a 2D z-buffer for depth testing. The
        z-buffer is crucial to ensure that only the parts of the torus closest to the observer are drawn.
    7.  Iterate through a grid of points on the surface of the original torus. The initial orientation is
        defined by parametric equations that place the torus hole along the Y-axis.
    8.  For each point:
        a. Apply the combined rotation matrix to transform the point to its new position.
        b. Project the 3D rotated point onto the 2D screen using an orthographic projection. This involves
           scaling and translating the x and y coordinates to fit the screen grid.
        c. Perform a depth check: if the z-coordinate of the new point is less than the value currently
           stored in the z-buffer at that screen location, it means this part of the surface is visible.
        d. If visible, update the z-buffer with the new depth and determine the appropriate character
           ('░', '▒', '▓', '█') based on the depth. Lighter shades are closer, darker shades are farther.
        e. Place the selected character into the character buffer.
    9.  After processing all points, print the resulting character buffer to the console, forming the ASCII image.
    10. Finally, print the original rotation values as requested.
    """
    # Parameters for the torus and the screen
    R1 = 2.0  # Major radius
    R2 = 1.0  # Minor radius
    screen_width = 50
    screen_height = 25

    # Rotation angles in degrees
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35

    # Shading characters from nearest to farthest
    chars = "░▒▓█"

    # Convert degrees to radians and negate for clockwise rotation
    rot_x_rad = -math.radians(rot_x_deg)
    rot_y_rad = -math.radians(rot_y_deg)
    rot_z_rad = -math.radians(rot_z_deg)

    # Precompute sines and cosines for rotation matrices
    cosA, sinA = math.cos(rot_x_rad), math.sin(rot_x_rad)
    cosB, sinB = math.cos(rot_y_rad), math.sin(rot_y_rad)
    cosC, sinC = math.cos(rot_z_rad), math.sin(rot_z_rad)

    # Rotation matrices
    Rx = np.array([
        [1, 0, 0],
        [0, cosA, -sinA],
        [0, sinA, cosA]
    ])
    Ry = np.array([
        [cosB, 0, sinB],
        [0, 1, 0],
        [-sinB, 0, cosB]
    ])
    Rz = np.array([
        [cosC, -sinC, 0],
        [sinC, cosC, 0],
        [0, 0, 1]
    ])

    # Combined rotation matrix (order: X, then Y, then Z)
    R_combined = Rz @ Ry @ Rx

    # Initialize buffers
    output_buffer = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    z_buffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # Define the density of points to draw on the torus surface
    theta_spacing = 0.07
    phi_spacing = 0.02

    # Generate and project torus points
    for theta in np.arange(0, 2 * math.pi, theta_spacing):
        for phi in np.arange(0, 2 * math.pi, phi_spacing):
            # Parametric equations for a torus with its hole along the Y-axis
            costheta, sintheta = math.cos(theta), math.sin(theta)
            cosphi, sinphi = math.cos(phi), math.sin(phi)
            
            circlex = R1 + R2 * costheta
            
            x = circlex * cosphi
            y = R2 * sintheta
            z = circlex * sinphi
            
            initial_point = np.array([x, y, z])
            
            # Rotate the point
            rotated_point = R_combined @ initial_point
            x_rot, y_rot, z_rot = rotated_point[0], rotated_point[1], rotated_point[2]

            # Project the 3D point to 2D screen coordinates
            # Scale factor to fit the torus within the screen
            scale = 10
            # Adjust y-scale for character aspect ratio (characters are taller than they are wide)
            y_scale_adjust = 0.5 
            x_proj = int(screen_width / 2 + scale * x_rot)
            y_proj = int(screen_height / 2 - scale * y_rot * y_scale_adjust)

            # Check if the projected point is within screen bounds and is the closest so far
            if 0 <= y_proj < screen_height and 0 <= x_proj < screen_width:
                if z_rot < z_buffer[y_proj][x_proj]:
                    z_buffer[y_proj][x_proj] = z_rot
                    
                    # Determine shading character based on depth (z_rot)
                    # Normalize z from ~-(R1+R2) to +(R1+R2) into an index for the chars string
                    total_radius = R1 + R2
                    luminance_index = int(((z_rot + total_radius) / (2 * total_radius)) * len(chars))
                    luminance_index = max(0, min(luminance_index, len(chars) - 1))
                    
                    output_buffer[y_proj][x_proj] = chars[luminance_index]

    # Print the final rendered image
    for row in output_buffer:
        print("".join(row))

    # Print the rotation values used in the calculation
    print(f"Rotation(X, Y, Z) = ({rot_x_deg}, {rot_y_deg}, {rot_z_deg})")

# Execute the function to generate the output
solve_torus_rotation()