import numpy as np

def solve():
    """
    Generates an ASCII representation of a rotated torus.
    This function models a torus, applies 3D rotation, calculates lighting,
    and renders the result to the console using a z-buffer.
    """

    # --- Model and Rendering Parameters ---
    R, r = 2.0, 0.9  # Major and minor radii of the torus
    width, height = 50, 22 # Screen dimensions
    screen_scale = 8 # Scale factor for projection
    light_source = np.array([0, 0, -1])  # Light from observer's direction
    shades = ['░', '▒', '▓', '█']

    # --- Rotation Angles ---
    # The problem specifies rotation in (X, Y, Z) degrees.
    rot_x_deg = 140
    rot_y_deg = 75
    rot_z_deg = 35

    # Print the equation as requested
    print(f"Equation for the rotation: ({rot_x_deg}, {rot_y_deg}, {rot_z_deg})")
    
    # --- Setup ---
    # Initialize screen grid and depth buffer
    screen = [[' ' for _ in range(width)] for _ in range(height)]
    z_buffer = [[float('inf') for _ in range(width)] for _ in range(height)]

    # Convert degrees to radians and negate for clockwise rotation
    ax = np.radians(-rot_x_deg)
    ay = np.radians(-rot_y_deg)
    az = np.radians(-rot_z_deg)

    # Create rotation matrices
    Rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
    Ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
    Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])

    # Combined rotation matrix (apply in order X, then Y, then Z)
    rotation_M = Rz @ Ry @ Rx

    # --- Main Loop ---
    # Iterate over the torus surface using parametric angles u and v
    u_steps = 120
    v_steps = 60
    for u in np.linspace(0, 2 * np.pi, u_steps):
        for v in np.linspace(0, 2 * np.pi, v_steps):
            # Original point P on a torus standing in the XZ-plane
            x = (R + r * np.cos(v)) * np.cos(u)
            y = r * np.sin(v)
            z = (R + r * np.cos(v)) * np.sin(u)
            P = np.array([x, y, z])

            # Original normal vector N for lighting calculation
            Nx = np.cos(v) * np.cos(u)
            Ny = np.sin(v)
            Nz = np.cos(v) * np.sin(u)
            N = np.array([Nx, Ny, Nz])

            # Apply rotation to the point and its normal vector
            P_rotated = rotation_M @ P
            N_rotated = rotation_M @ N

            # Calculate illumination (dot product of normal and light source)
            illumination = np.dot(N_rotated, light_source)

            # Render only the surfaces facing the light source
            if illumination > 0:
                px, py, pz = P_rotated

                # Project the 3D point to 2D screen coordinates
                screen_x = int(width / 2 + px * screen_scale)
                # Invert y-coordinate for screen (row 0 is top)
                screen_y = int(height / 2 - py * screen_scale)

                # Z-buffer check for occlusion
                if 0 <= screen_x < width and 0 <= screen_y < height:
                    if pz < z_buffer[screen_y][screen_x]:
                        z_buffer[screen_y][screen_x] = pz
                        
                        # Determine which character to use based on illumination
                        shade_index = min(int(illumination * len(shades)), len(shades) - 1)
                        screen[screen_y][screen_x] = shades[shade_index]

    # --- Print the final rendered image ---
    for row in screen:
        print("".join(row))

solve()