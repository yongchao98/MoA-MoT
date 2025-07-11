import math

def generate_rotated_torus_view():
    """
    This function simulates and renders a rotated torus based on the problem description.
    """
    
    # 1. Define parameters
    R1 = 2  # Major radius
    R2 = 1  # Minor radius
    width, height = 55, 25  # Screen dimensions

    # Rotation angles from the problem, converted to radians.
    # The problem states positive rotation is clockwise, so we negate the angles
    # for use in standard counter-clockwise rotation matrices.
    rot_x_deg, rot_y_deg, rot_z_deg = 140, 75, 35
    
    # Initial rotation to get the 'standing up' torus view.
    # Based on analysis, this must be a standard CCW rotation to match the initial image shading.
    initial_rot_y_deg = 90 

    ax = math.radians(-rot_x_deg)
    ay = math.radians(-rot_y_deg)
    az = math.radians(-rot_z_deg)
    initial_ay = math.radians(initial_rot_y_deg)

    print("Thinking Process & Plan:")
    print("1. Model a torus in 3D space with its axis of revolution on the Z-axis.")
    print("2. Apply an initial 90-degree rotation around the Y-axis to match the provided starting view.")
    print("3. Sequentially apply the specified new rotations around the X, Y, and Z axes.")
    print(f"   - Rotation around X: {rot_x_deg} degrees (clockwise)")
    print(f"   - Rotation around Y: {rot_y_deg} degrees (clockwise)")
    print(f"   - Rotation around Z: {rot_z_deg} degrees (clockwise)")
    print("4. Project the final 3D points onto a 2D screen using a z-buffer for occlusion.")
    print("5. Shade the resulting view based on depth, where lighter shades are nearer.")
    print("\nCalculating final view...\n")


    # 2. Matrix multiplication helpers
    def multiply_mat_vec(matrix, vector):
        res = [0, 0, 0]
        for i in range(3):
            for j in range(3):
                res[i] += matrix[i][j] * vector[j]
        return res

    def multiply_mat_mat(m1, m2):
        res = [[0] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    res[i][j] += m1[i][k] * m2[k][j]
        return res

    # 3. Define and combine rotation matrices
    Rx = [
        [1, 0, 0],
        [0, math.cos(ax), -math.sin(ax)],
        [0, math.sin(ax), math.cos(ax)]
    ]
    Ry = [
        [math.cos(ay), 0, math.sin(ay)],
        [0, 1, 0],
        [-math.sin(ay), 0, math.cos(ay)]
    ]
    Rz = [
        [math.cos(az), -math.sin(az), 0],
        [math.sin(az), math.cos(az), 0],
        [0, 0, 1]
    ]
    R_initial = [
        [math.cos(initial_ay), 0, math.sin(initial_ay)],
        [0, 1, 0],
        [-math.sin(initial_ay), 0, math.cos(initial_ay)]
    ]

    # Total transformation: Rz @ Ry @ Rx @ R_initial
    R_new = multiply_mat_mat(Rz, multiply_mat_mat(Ry, Rx))
    R_total = multiply_mat_mat(R_new, R_initial)

    # 4. Render the torus
    output = [[' ' for _ in range(width)] for _ in range(height)]
    zbuffer = [[float('inf') for _ in range(width)] for _ in range(height)]
    
    chars = "░▒▓█" # Nearer (light) to farther (dark)
    theta_spacing = 0.07
    phi_spacing = 0.025
    
    # Loop through torus surface points
    theta = 0
    while theta < 2 * math.pi:
        phi = 0
        while phi < 2 * math.pi:
            # Base torus point (axis on Z)
            cos_phi, sin_phi = math.cos(phi), math.sin(phi)
            cos_theta, sin_theta = math.cos(theta), math.sin(theta)
            
            x = (R1 + R2 * cos_theta) * cos_phi
            y = (R1 + R2 * cos_theta) * sin_phi
            z = R2 * sin_theta
            
            # Apply total rotation
            p_rotated = multiply_mat_vec(R_total, [x, y, z])
            xr, yr, zr = p_rotated[0], p_rotated[1], p_rotated[2]

            # Project to screen
            screen_scale = width / (2.5 * (R1 + R2)) * 0.8
            # Compensate for character aspect ratio
            xp = int(width / 2 + screen_scale * xr * 2) 
            yp = int(height / 2 - screen_scale * yr)

            if 0 <= xp < width and 0 <= yp < height:
                if zr < zbuffer[yp][xp]:
                    zbuffer[yp][xp] = zr
            phi += phi_spacing
        theta += theta_spacing
    
    # 5. Determine shading from z-buffer
    min_z = min(z for row in zbuffer for z in row if z != float('inf'))
    max_z = max(z for row in zbuffer for z in row if z != float('inf'))
    z_range = max_z - min_z
    if z_range == 0: z_range = 1 # Avoid division by zero

    for r in range(height):
        for c in range(width):
            if zbuffer[r][c] != float('inf'):
                z = zbuffer[r][c]
                normalized_z = (z - min_z) / z_range
                char_index = int(normalized_z * (len(chars) - 1))
                output[r][c] = chars[char_index]
    
    # 6. Print the result
    for row in output:
        print("".join(row))

generate_rotated_torus_view()
<<<C>>>