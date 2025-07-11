import math

def solve_torus_rotation():
    """
    This script simulates the rotation of a 3D torus and renders it as ASCII art.
    """
    # --- Parameters ---
    R, r = 2.0, 1.0  # Major and minor radii of the torus
    screen_width, screen_height = 52, 22 # Canvas dimensions to match options
    
    # Scaling factors for projection
    x_scale = 7
    y_scale = 14 # Account for character aspect ratio
    
    # Angle steps for drawing the torus surface
    theta_step, phi_step = 0.04, 0.04

    # --- Rotation Angles ---
    rx_deg, ry_deg, rz_deg = 140, 75, 35
    rx = math.radians(rx_deg)
    ry = math.radians(ry_deg)
    rz = math.radians(rz_deg)
    
    print(f"Applying clockwise rotations (X, Y, Z): ({rx_deg}, {ry_deg}, {rz_deg}) degrees\n")


    # --- Clockwise Rotation Matrices ---
    cos_rx, sin_rx = math.cos(rx), math.sin(rx)
    Rx_cw = [[1, 0, 0], [0, cos_rx, sin_rx], [0, -sin_rx, cos_rx]]

    cos_ry, sin_ry = math.cos(ry), math.sin(ry)
    Ry_cw = [[cos_ry, 0, -sin_ry], [0, 1, 0], [sin_ry, 0, cos_ry]]

    cos_rz, sin_rz = math.cos(rz), math.sin(rz)
    Rz_cw = [[cos_rz, sin_rz, 0], [-sin_rz, cos_rz, 0], [0, 0, 1]]

    # --- Matrix Multiplication Helper ---
    def matmul(A, B):
        return [[sum(a * b for a, b in zip(A_row, B_col)) for B_col in zip(*B)] for A_row in A]

    # Combined rotation matrix: R_final = Rz * Ry * Rx
    R_total = matmul(Rz_cw, matmul(Ry_cw, Rx_cw))

    # --- Buffers ---
    canvas = [[' ' for _ in range(screen_width)] for _ in range(screen_height)]
    zbuffer = [[float('inf') for _ in range(screen_width)] for _ in range(screen_height)]

    # --- Torus Generation and Rendering Loop ---
    for theta in range(int(2 * math.pi / theta_step)):
        theta_val = theta * theta_step
        cos_theta, sin_theta = math.cos(theta_val), math.sin(theta_val)
        
        for phi in range(int(2 * math.pi / phi_step)):
            phi_val = phi * phi_step
            cos_phi, sin_phi = math.cos(phi_val), math.sin(phi_val)

            # Initial point on a torus with hole along the y-axis
            x0 = (R + r * cos_phi) * cos_theta
            y0 = r * sin_phi
            z0 = -(R + r * cos_phi) * sin_theta
            p0 = [x0, y0, z0]
            
            # Initial normal vector
            nx0 = cos_phi * cos_theta
            ny0 = sin_phi
            nz0 = -cos_phi * sin_theta
            n0 = [nx0, ny0, nz0]

            # Apply rotation by multiplying with the combined matrix
            p_rot = [sum(R_total[i][j] * p0[j] for j in range(3)) for i in range(3)]
            n_rot = [sum(R_total[i][j] * n0[j] for j in range(3)) for i in range(3)]
            
            # Back-face culling (only render if normal's z-component faces camera)
            if n_rot[2] < 0:
                x, y, z = p_rot[0], p_rot[1], p_rot[2]
                
                # Orthographic projection onto screen
                sx = int(screen_width / 2 + x * x_scale)
                sy = int(screen_height / 2 - y * y_scale / 2.0)

                if 0 <= sx < screen_width and 0 <= sy < screen_height:
                    if z < zbuffer[sy][sx]:
                        zbuffer[sy][sx] = z

    # --- Find Z-range for shading and fill canvas ---
    min_z, max_z = float('inf'), float('-inf')
    for r_idx in range(screen_height):
        for c_idx in range(screen_width):
            if zbuffer[r_idx][c_idx] != float('inf'):
                min_z = min(min_z, zbuffer[r_idx][c_idx])
                max_z = max(max_z, zbuffer[r_idx][c_idx])

    shades = ['░', '▒', '▓', '█']
    z_range = max_z - min_z if (max_z - min_z) != 0 else 1

    for r_idx in range(screen_height):
        for c_idx in range(screen_width):
            z = zbuffer[r_idx][c_idx]
            if z != float('inf'):
                shade_index = int(((z - min_z) / z_range) * 3.99)
                canvas[r_idx][c_idx] = shades[shade_index]

    # --- Print final result ---
    for row in canvas:
        print("".join(row))

# Execute the simulation
solve_torus_rotation()