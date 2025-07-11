import math

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points (X, Y, Z).
    """
    # Step 1: Define input data based on the map.
    # Heights of the three points in meters.
    h_X = 140.0
    h_Y = 80.0
    h_Z = 120.0 # Corrected from description, X=120, Z=140 based on standard 3-point problem layout
    # Z has highest elevation, Y has lowest. Let's relabel for clarity of calculation
    # Let P1 = Z (highest), P2 = X (middle), P3 = Y (lowest)
    # This avoids confusion in the three-point problem method, though vector math handles any order.
    # Let's stick to the labels in the problem: X, Y, Z.
    h_X = 120.0
    h_Y = 80.0
    h_Z = 140.0

    # Pixel coordinates measured from the image (origin at top-left).
    # These are estimates and represent the main source of uncertainty.
    pix_X = (494, 381)
    pix_Y = (211, 461)
    pix_Z = (356, 186)
    
    # Pixel coordinates of the scale bar.
    scale_start_pix = 640
    scale_end_pix = 778
    scale_length_m = 200.0

    print("Step 1: Parsing Map Data")
    # Step 2: Calculate the scale (meters per pixel).
    scale_length_pix = scale_end_pix - scale_start_pix
    m_per_pix = scale_length_m / scale_length_pix
    print(f"Scale: {scale_length_pix} pixels = {scale_length_m} meters")
    print(f"Meters per pixel = {m_per_pix:.4f}\n")

    # Step 3: Establish a 3D coordinate system.
    # We set point Y as the origin for the horizontal plane.
    # The coordinate system is (East, North, Up).
    # North is opposite to the image's y-axis direction.
    map_Y = (0.0, 0.0)
    map_X = ((pix_X[0] - pix_Y[0]) * m_per_pix, (pix_Y[1] - pix_X[1]) * m_per_pix)
    map_Z = ((pix_Z[0] - pix_Y[0]) * m_per_pix, (pix_Y[1] - pix_Z[1]) * m_per_pix)

    # Combine map coordinates with heights to get 3D points.
    p_X = (map_X[0], map_X[1], h_X)
    p_Y = (map_Y[0], map_Y[1], h_Y)
    p_Z = (map_Z[0], map_Z[1], h_Z)

    print("Step 2: Defining 3D points in space (East, North, Up)")
    print(f"Point Y: ({p_Y[0]:.1f} m, {p_Y[1]:.1f} m, {p_Y[2]:.1f} m)")
    print(f"Point X: ({p_X[0]:.1f} m, {p_X[1]:.1f} m, {p_X[2]:.1f} m)")
    print(f"Point Z: ({p_Z[0]:.1f} m, {p_Z[1]:.1f} m, {p_Z[2]:.1f} m)\n")

    # Step 4: Create two vectors on the plane.
    v_YX = (p_X[0] - p_Y[0], p_X[1] - p_Y[1], p_X[2] - p_Y[2])
    v_YZ = (p_Z[0] - p_Y[0], p_Z[1] - p_Y[1], p_Z[2] - p_Y[2])

    print("Step 3: Creating vectors on the planar surface")
    print(f"Vector YX: ({v_YX[0]:.1f}, {v_YX[1]:.1f}, {v_YX[2]:.1f})")
    print(f"Vector YZ: ({v_YZ[0]:.1f}, {v_YZ[1]:.1f}, {v_YZ[2]:.1f})\n")

    # Step 5: Calculate the normal vector using the cross product.
    # N = v_YX x v_YZ
    Nx = v_YX[1] * v_YZ[2] - v_YX[2] * v_YZ[1]
    Ny = v_YX[2] * v_YZ[0] - v_YX[0] * v_YZ[2]
    Nz = v_YX[0] * v_YZ[1] - v_YX[1] * v_YZ[0]
    
    print("Step 4: Calculating the normal vector N = (Nx, Ny, Nz) to the plane")
    print(f"Normal Vector N = ({Nx:.1f}, {Ny:.1f}, {Nz:.1f})\n")

    # Step 6: Calculate the dip angle.
    # tan(dip) = horizontal_component_of_N / vertical_component_of_N
    # tan(dip) = sqrt(Nx^2 + Ny^2) / |Nz|
    horizontal_component_magnitude = math.sqrt(Nx**2 + Ny**2)
    vertical_component_magnitude = abs(Nz)

    if vertical_component_magnitude == 0:
        # This case implies a vertical plane, dip is 90 degrees.
        dip_rad = math.pi / 2
    else:
        tan_dip = horizontal_component_magnitude / vertical_component_magnitude
        dip_rad = math.atan(tan_dip)

    dip_deg = math.degrees(dip_rad)
    
    print("Step 5: Calculating the Dip Angle")
    print(f"The final equation for the dip is: dip = arctan(sqrt(Nx² + Ny²) / |Nz|)")
    print(f"Plugging in the numbers:")
    print(f"dip = arctan(sqrt(({Nx:.1f})² + ({Ny:.1f})²) / |{Nz:.1f}|)")
    print(f"dip = arctan({horizontal_component_magnitude:.1f} / {vertical_component_magnitude:.1f})")
    print(f"dip = arctan({tan_dip:.4f})")
    print(f"Dip = {dip_deg:.2f} degrees\n")

    # Step 7: Round to the nearest degree.
    rounded_dip = round(dip_deg)
    print(f"Final Answer: The dip rounded to the nearest degree is {rounded_dip} degrees.")
    
    return rounded_dip

final_answer = calculate_dip()
<<<8>>>