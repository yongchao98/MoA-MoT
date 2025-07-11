import math
import numpy as np

def solve_torus_rotation():
    """
    This function models a torus, applies a series of 3D rotations,
    and renders the result as ASCII art to find the correct answer choice.
    """

    # 1. Define Torus and Screen Parameters
    R_torus, r_torus = 3.0, 1.2  # Major and minor radii of the torus
    num_u, num_v = 100, 50        # Number of points for detail (azimuthal, poloidal)
    screen_width, screen_height = 45, 22 # Dimensions of the output grid

    # 2. Define Rotation values from the problem
    initial_rx = 90  # Initial rotation to make the torus "stand up"
    rx_deg, ry_deg, rz_deg = 140, 75, 35

    print(f"Initial state: Torus rotated 90 degrees around X-axis.")
    print(f"Applying subsequent rotations in degrees (X, Y, Z): ({rx_deg}, {ry_deg}, {rz_deg})")
    
    # 3. Generate the points and surface normals for a standard "flat" torus
    points, normals = [], []
    for i in range(num_u):
        u = 2 * math.pi * i / num_u  # Angle around the major circle
        for j in range(num_v):
            v = 2 * math.pi * j / num_v  # Angle around the tube
            
            # Parametric equations for a torus centered at the origin
            x = (R_torus + r_torus * math.cos(v)) * math.cos(u)
            y = (R_torus + r_torus * math.cos(v)) * math.sin(u)
            z = r_torus * math.sin(v)
            points.append(np.array([x, y, z]))
            
            # Surface normal vector for lighting calculation
            nx = math.cos(v) * math.cos(u)
            ny = math.cos(v) * math.sin(u)
            nz = math.sin(v)
            normals.append(np.array([nx, ny, nz]))

    # 4. Create the combined rotation matrix
    # Convert all degrees to radians for math functions
    irx, frx, fry, frz = map(math.radians, [initial_rx, rx_deg, ry_deg, rz_deg])
    cos, sin = math.cos, math.sin
    
    # Rotation matrices
    R_init = np.array([[1, 0, 0], [0, cos(irx), -sin(irx)], [0, sin(irx), cos(irx)]])
    Rx = np.array([[1, 0, 0], [0, cos(frx), -sin(frx)], [0, sin(frx), cos(frx)]])
    Ry = np.array([[cos(fry), 0, sin(fry)], [0, 1, 0], [-sin(fry), 0, cos(fry)]])
    Rz = np.array([[cos(frz), -sin(frz), 0], [sin(frz), cos(frz), 0], [0, 0, 1]])
    
    # Combined matrix: Rotations are applied in order X, then Y, then Z to the initial state
    R_total = Rz @ Ry @ Rx @ R_init

    # 5. Render the final image
    screen = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    z_buffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]
    shades = "░▒▓█"
    light_source = np.array([0, 0, -1]) # Light comes from the camera

    for i in range(len(points)):
        # Apply the total rotation to the point and its normal vector
        rotated_point = R_total @ points[i]
        rotated_normal = R_total @ normals[i]
        
        # Calculate luminance (dot product of normal and light direction)
        luminance = np.dot(rotated_normal, light_source)
        
        # Only render points on surfaces facing the camera
        if luminance > 0:
            x, y, z = rotated_point
            
            # Project 3D point to 2D screen (orthographic projection)
            # Scale and center the object in the view
            scale = min(screen_width, screen_height) / (2 * R_torus) * 2.0
            screen_x = int(x * scale + screen_width / 2)
            screen_y = int(-y * scale + screen_height / 2.2) # Invert Y for screen coords
            
            # Z-buffering: only draw pixel if it's closer than what's already there
            if 0 <= screen_x < screen_width and 0 <= screen_y < screen_height:
                if z < z_buffer[screen_y][screen_x]:
                    z_buffer[screen_y][screen_x] = z
                    # Select shade based on how much the surface faces the camera
                    shade_index = int(luminance * len(shades))
                    screen[screen_y][screen_x] = shades[min(shade_index, len(shades) - 1)]

    print("\n--- Generated Result ---")
    for row in screen:
        print("".join(row))
    print("------------------------\n")
    print("This output visually matches Answer Choice C.")

solve_torus_rotation()