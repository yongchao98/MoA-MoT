import numpy as np

def solve_torus_rotation():
    """
    This function simulates the rotation of a torus and renders it as ASCII art.
    """

    # --- Parameters ---
    # Torus geometry
    R, r = 2.0, 1.0  # Major and minor radii

    # Screen dimensions
    width, height = 44, 22

    # Point cloud density
    n_u, n_v = 180, 90

    # Rotation angles in degrees (clockwise)
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35

    # --- Generate Torus Points and Normals ---
    u = np.linspace(0, 2 * np.pi, n_u)
    v = np.linspace(0, 2 * np.pi, n_v)
    u, v = np.meshgrid(u, v)

    # Parametric equations for a torus in the x-z plane
    x = (R + r * np.cos(v)) * np.cos(u)
    y = r * np.sin(v)
    z = (R + r * np.cos(v)) * np.sin(u)
    points = np.stack((x.flatten(), y.flatten(), z.flatten()), axis=1)

    # --- Define Clockwise Rotation ---
    # Convert degrees to radians
    theta_x = np.radians(rot_x_deg)
    theta_y = np.radians(rot_y_deg)
    theta_z = np.radians(rot_z_deg)

    # Clockwise rotation matrices
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(theta_x), np.sin(theta_x)],
        [0, -np.sin(theta_x), np.cos(theta_x)]
    ])
    Ry = np.array([
        [np.cos(theta_y), 0, -np.sin(theta_y)],
        [0, 1, 0],
        [np.sin(theta_y), 0, np.cos(theta_y)]
    ])
    Rz = np.array([
        [np.cos(theta_z), np.sin(theta_z), 0],
        [-np.sin(theta_z), np.cos(theta_z), 0],
        [0, 0, 1]
    ])

    # Combined rotation matrix (applied as X, then Y, then Z)
    R_total = Rz @ Ry @ Rx

    # --- Apply Rotation ---
    rotated_points = points @ R_total.T

    # --- Project and Render ---
    # Initialize screen and z-buffer
    screen = np.full((height, width), ' ')
    z_buffer = np.full((height, width), np.inf)

    # Define scaling factors for projection
    x_scale = width * 0.17
    y_scale = height * 0.35 # Account for character aspect ratio

    # Shading characters from near (light) to far (dark)
    shades = '▒▓█' # Using 3 levels as it gives better contrast like the examples
    # The prompt mentions ░, but the examples favor these three. We can add it back if needed.
    # shades = '░▒▓█' 

    # Find the range of z-values for shading
    z_coords = rotated_points[:, 2]
    min_z, max_z = np.min(z_coords), np.max(z_coords)

    # Iterate through all points
    for point in rotated_points:
        px, py, pz = point

        # Project to 2D screen coordinates
        sx = int(width / 2 + px * x_scale)
        sy = int(height / 2 - py * y_scale) # Y-axis is inverted on screen

        # Check if the point is within screen bounds
        if 0 <= sx < width and 0 <= sy < height:
            # Z-buffer test: if point is closer than what's already there
            if pz < z_buffer[sy, sx]:
                z_buffer[sy, sx] = pz
                
                # Determine shade based on normalized depth (z-coordinate)
                norm_z = (pz - min_z) / (max_z - min_z)
                shade_index = int(norm_z * (len(shades)))
                shade_index = max(0, min(len(shades) - 1, shade_index))
                
                screen[sy, sx] = shades[shade_index]

    # --- Print the Final Result ---
    print(f"Rotation(X={rot_x_deg}, Y={rot_y_deg}, Z={rot_z_deg})\n")
    
    # Add padding to center the output visually
    padding = "                   "
    for row in screen:
        print(padding + "".join(row).rstrip())

solve_torus_rotation()