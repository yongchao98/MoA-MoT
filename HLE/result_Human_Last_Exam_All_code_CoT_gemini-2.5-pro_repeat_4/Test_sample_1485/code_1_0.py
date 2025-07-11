import math

def generate_rotated_torus_view():
    """
    This script simulates the 3D rotation of a torus and projects it onto a 2D screen.
    The final ASCII art is generated based on the depth of the torus surface.
    """
    # 1. Define torus shape and screen parameters.
    #    R1 is the major radius (center of torus to center of tube).
    #    R2 is the minor radius (radius of the tube).
    R1 = 2.0
    R2 = 1.0
    screen_width = 45
    screen_height = 25
    screen_scale = 10 # Controls the size of the torus on screen

    # 2. Define the rotation angles from the problem description.
    #    The problem implies a clockwise positive rotation.
    angle_x_deg, angle_y_deg, angle_z_deg = 140, 75, 35

    # 3. Initialize buffers.
    #    - screen: A 2D array to hold the final characters.
    #    - zbuffer: A 2D array to handle occlusion (depth). Closer points overwrite farther ones.
    screen = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # 4. Convert angles to radians and pre-calculate sines and cosines for efficiency.
    ax = math.radians(angle_x_deg)
    ay = math.radians(angle_y_deg)
    az = math.radians(angle_z_deg)
    cosA, sinA = math.cos(ax), math.sin(ax)
    cosB, sinB = math.cos(ay), math.sin(ay)
    cosC, sinC = math.cos(az), math.sin(az)

    # 5. First Pass: Iterate through the surface of the torus.
    #    For each point, apply the rotations and project it to the screen.
    #    Store the depth (z-coordinate) in the z-buffer.
    theta_step = 0.07  # Granularity of the tube's circle
    phi_step = 0.02    # Granularity of the torus's main circle

    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            # Generate the initial point on a torus lying in the XZ plane (hole along Y-axis).
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            circlex = R2 * cos_theta + R1
            x = circlex * cos_phi
            y = R2 * sin_theta
            z = circlex * sin_phi

            # Apply rotations sequentially (X, then Y, then Z) for clockwise rotation.
            # Rotate around X-axis
            y1 = y * cosA + z * sinA
            z1 = -y * sinA + z * cosA
            x1 = x
            
            # Rotate around Y-axis
            x2 = x1 * cosB - z1 * sinB
            z2 = x1 * sinB + z1 * cosB
            y2 = y1

            # Rotate around Z-axis
            x3 = x2 * cosC + y2 * sinC
            y3 = -x2 * sinC + y2 * cosC
            z3 = z2
            
            # Project the 3D point to 2D screen coordinates.
            # The screen's Y-axis is inverted (0 is at the top).
            xp = int(screen_width / 2 + screen_scale * x3)
            yp = int(screen_height / 2 - screen_scale * y3)

            # Update the z-buffer if this point is closer to the observer.
            if 0 <= xp < screen_width and 0 <= yp < screen_height:
                if z3 < zbuffer[yp][xp]:
                    zbuffer[yp][xp] = z3
            
            phi += phi_step
        theta += theta_step

    # 6. Second Pass: Determine shading from z-buffer and draw.
    #    Find the min and max depth of the visible points to normalize depth values.
    min_z = float('inf')
    max_z = float('-inf')
    for r in zbuffer:
        for val in r:
            if val != float('inf'):
                min_z = min(min_z, val)
                max_z = max(max_z, val)

    z_range = max_z - min_z if max_z > min_z else 1.0

    # The shades correspond to depth, from nearest (░) to farthest (█).
    shades = ['░', '▒', '▓', '█']

    for r in range(screen_height):
        for c in range(screen_width):
            z_val = zbuffer[r][c]
            if z_val != float('inf'):
                # Normalize the depth to a value between 0 and 1.
                normalized_z = (z_val - min_z) / z_range
                # Map the normalized depth to a shade character.
                shade_index = int(normalized_z * 3.99)
                screen[r][c] = shades[shade_index]

    # 7. Print the final result.
    # Output the rotation values used as per the instructions.
    print(f"Rotation(X, Y, Z) = ({angle_x_deg}, {angle_y_deg}, {angle_z_deg})")
    for row in screen:
        print("".join(row))

generate_rotated_torus_view()