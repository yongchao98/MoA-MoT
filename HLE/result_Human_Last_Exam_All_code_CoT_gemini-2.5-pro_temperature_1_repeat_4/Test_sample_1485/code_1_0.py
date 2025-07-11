import numpy as np
import math

def generate_rotated_torus_view():
    """
    Generates an ASCII representation of a torus after specified rotations.
    """
    # Rotation angles in degrees (X, Y, Z)
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35

    # Convert degrees to radians for math functions.
    # Positive rotation is clockwise, so we use negative angles with standard CCW rotation matrices.
    theta_x = math.radians(-rot_x_deg)
    theta_y = math.radians(-rot_y_deg)
    theta_z = math.radians(-rot_z_deg)

    # Torus parameters
    R, r = 2.0, 1.0  # Major and minor radii

    # Screen parameters
    width, height = 50, 22
    
    # Create rotation matrices
    Rx = np.array([[1, 0, 0], [0, math.cos(theta_x), -math.sin(theta_x)], [0, math.sin(theta_x), math.cos(theta_x)]])
    Ry = np.array([[math.cos(theta_y), 0, math.sin(theta_y)], [0, 1, 0], [-math.sin(theta_y), 0, math.cos(theta_y)]])
    Rz = np.array([[math.cos(theta_z), -math.sin(theta_z), 0], [math.sin(theta_z), math.cos(theta_z), 0], [0, 0, 1]])
    
    # Combined rotation matrix (applied in order X, then Y, then Z)
    R_combined = Rz @ Ry @ Rx

    # Buffers for rendering
    z_buffer = [[-float('inf')] * width for _ in range(height)]
    # Store the true z-value for shading
    z_values_for_shading = [[None] * width for _ in range(height)]

    # First Pass: Generate points, rotate, project, and fill z-buffer
    u_step, v_step = 0.04, 0.04
    u = 0
    while u < 2 * math.pi:
        v = 0
        while v < 2 * math.pi:
            # Parametrization for a torus initially revolving around the Y-axis
            x0 = (R + r * math.cos(u)) * math.sin(v)
            y0 = r * math.sin(u)
            z0 = (R + r * math.cos(u)) * math.cos(v)
            P_initial = np.array([x0, y0, z0])
            
            # Normal vector for back-face culling
            nx0 = math.cos(u) * math.sin(v)
            ny0 = math.sin(u)
            nz0 = math.cos(u) * math.cos(v)
            N_initial = np.array([nx0, ny0, nz0])

            # Apply the combined rotation
            P_rotated = R_combined @ P_initial
            N_rotated = R_combined @ N_initial
            
            # Check if the surface is facing the observer (view is along -Z axis)
            if np.dot(N_rotated, np.array([0, 0, -1])) > 0:
                # Add a distance offset to push the object away from the camera at (0,0,0)
                z_final = P_rotated[2] + 5.0
                
                # Project onto 2D screen (perspective projection)
                ooz = 1.0 / z_final if z_final != 0 else 0
                scale = 22
                xp = int(width / 2 + scale * P_rotated[0] * ooz)
                yp = int(height / 2 - scale * P_rotated[1] * ooz)
                
                if 0 <= xp < width and 0 <= yp < height:
                    # Z-buffer check: is this point closer than the one already there?
                    if ooz > z_buffer[yp][xp]:
                        z_buffer[yp][xp] = ooz
                        # Store the point's original z-coordinate for shading
                        z_values_for_shading[yp][xp] = P_rotated[2]
            v += v_step
        u += u_step

    # Second Pass: Determine min/max Z for shading and render the final screen
    min_z, max_z = float('inf'), -float('inf')
    for row in z_values_for_shading:
        for z in row:
            if z is not None:
                min_z = min(min_z, z)
                max_z = max(max_z, z)

    z_range = max_z - min_z if max_z > min_z else 1.0
    shades = "░▒▓█"
    screen = [[' '] * width for _ in range(height)]

    for r_idx in range(height):
        for c_idx in range(width):
            if z_values_for_shading[r_idx][c_idx] is not None:
                # Normalize z from 0 to 1 (0=nearest, 1=farthest)
                norm_z = (z_values_for_shading[r_idx][c_idx] - min_z) / z_range
                # Map normalized z to shade character
                shade_index = min(int(norm_z * len(shades)), len(shades) - 1)
                screen[r_idx][c_idx] = shades[shade_index]

    # Print the final result
    print("The final view of the torus after rotation is calculated below.")
    print(f"The equation for the rotation is a sequential rotation around the X, Y, and Z axes by the angles ({rot_x_deg}, {rot_y_deg}, {rot_z_deg}) degrees.")
    print("\n" * 2) # Add some spacing
    for row in screen:
        print("".join(row).rstrip())
    print("\n" * 10) # Add some spacing

generate_rotated_torus_view()