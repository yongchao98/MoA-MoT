import math
import numpy as np

def solve_torus_rotation():
    """
    Calculates and prints the ASCII representation of a rotated torus.
    """
    # 1. Define rotation parameters as provided in the problem.
    given_rot_x_deg = 140
    given_rot_y_deg = 75
    given_rot_z_deg = 35

    print(f"Initial rotation requested (X, Y, Z degrees): ({given_rot_x_deg}, {given_rot_y_deg}, {given_rot_z_deg})")

    # 2. Setup simulation parameters
    width, height = 60, 22
    # Torus geometry
    R1 = 2.0  # Major radius
    R2 = 1.0  # Minor radius
    # Orthographic projection and display parameters
    scale = 14.5
    aspect_ratio_correction = 0.5 # To account for non-square characters in a terminal

    # 3. Define total rotation
    # The initial view is a torus standing up, which is a 90-degree rotation
    # around the X-axis from a standard torus lying flat in the X-Y plane.
    # We add the requested rotation to this initial orientation.
    total_rot_x_deg = 90.0 + given_rot_x_deg
    total_rot_y_deg = float(given_rot_y_deg)
    total_rot_z_deg = float(given_rot_z_deg)

    # Convert final angles to radians for trigonometric functions
    ax = math.radians(total_rot_x_deg)
    ay = math.radians(total_rot_y_deg)
    az = math.radians(total_rot_z_deg)

    # Create rotation matrices for each axis
    Rx = np.array([
        [1, 0, 0],
        [0, math.cos(ax), -math.sin(ax)],
        [0, math.sin(ax), math.cos(ax)]
    ])
    Ry = np.array([
        [math.cos(ay), 0, math.sin(ay)],
        [0, 1, 0],
        [-math.sin(ay), 0, math.cos(ay)]
    ])
    Rz = np.array([
        [math.cos(az), -math.sin(az), 0],
        [math.sin(az), math.cos(az), 0],
        [0, 0, 1]
    ])
    # Combine into a single transformation matrix. Order is Z, then Y, then X.
    M = np.dot(Rz, np.dot(Ry, Rx))

    # 4. Initialize buffers
    # The screen holds the final characters.
    screen = [[' '] * width for _ in range(height)]
    # The z-buffer stores the depth of the closest point for each pixel.
    # We use a right-handed system where +Z is towards the viewer, so we initialize to -infinity.
    zbuffer = [[float('-inf')] * width for _ in range(height)]

    # 5. Generate points on the torus, rotate, and project them
    theta_step = 0.07
    phi_step = 0.07
    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            # Parametric equations for a standard torus in the x-y plane
            costheta, sintheta = math.cos(theta), math.sin(theta)
            cosphi, sinphi = math.cos(phi), math.sin(phi)
            x0 = (R1 + R2 * costheta) * cosphi
            y0 = (R1 + R2 * costheta) * sinphi
            z0 = R2 * sintheta
            
            p0 = np.array([x0, y0, z0])
            
            # Apply the combined rotation
            p_rotated = np.dot(M, p0)
            x, y, z = p_rotated[0], p_rotated[1], p_rotated[2]
            
            # Project the 3D point to 2D screen coordinates (orthographic)
            xp = int(width / 2 + x * scale)
            yp = int(height / 2 - y * scale * aspect_ratio_correction)
            
            # If the point is on screen and closer than what's in the z-buffer, update
            if 0 <= xp < width and 0 <= yp < height:
                if z > zbuffer[yp][xp]:
                    zbuffer[yp][xp] = z
            
            phi += phi_step
        theta += theta_step

    # 6. Determine shading based on depth
    # Find the min and max depth values across all visible points
    min_z = float('inf')
    max_z = float('-inf')
    for r in range(height):
        for c in range(width):
            if zbuffer[r][c] != float('-inf'):
                min_z = min(min_z, zbuffer[r][c])
                max_z = max(max_z, zbuffer[r][c])
    
    # Shading characters: Nearer points (higher Z) get lighter characters.
    chars = ['█', '▓', '▒', '░']
    z_range = max_z - min_z
    if z_range == 0: z_range = 1

    # Fill the screen buffer with the correctly shaded characters
    for r in range(height):
        for c in range(width):
            if zbuffer[r][c] != float('-inf'):
                # Normalize depth from 0 (farthest) to 1 (nearest)
                norm_z = (zbuffer[r][c] - min_z) / z_range
                # Select character based on normalized depth
                char_index = int(norm_z * (len(chars) - 0.01))
                screen[r][c] = chars[char_index]

    # 7. Print the final rendered image
    print("\nFinal rendered view:")
    for row in screen:
        print("".join(row))

solve_torus_rotation()