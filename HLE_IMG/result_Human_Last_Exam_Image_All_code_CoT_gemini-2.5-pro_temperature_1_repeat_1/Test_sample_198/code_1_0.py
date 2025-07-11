import math

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points (X, Y, Z).
    The calculation is based on coordinates and heights derived from the provided map.
    """
    # Step 1: Define the data from the map and problem description.
    # Heights in meters
    h_X = 120
    h_Y = 80
    h_Z = 140

    # Pixel coordinates measured from the image (origin at top-left corner)
    # These are estimated from the image provided.
    px_X, py_X = (612, 558)
    px_Y, py_Y = (208, 658)
    px_Z, py_Z = (460, 258)

    # Pixel coordinates of the scale bar for '0 m' and '200 m' marks
    scale_0_px = 760
    scale_200_px = 957
    
    # Step 2: Determine the map scale.
    scale_dist_px = scale_200_px - scale_0_px
    scale_dist_m = 200.0  # meters
    m_per_pixel = scale_dist_m / scale_dist_px
    
    print(f"Map Scale Calculation:")
    print(f"  Pixel distance for 200m = {scale_dist_px} px")
    print(f"  Scale = {scale_dist_m} m / {scale_dist_px} px = {m_per_pixel:.4f} m/px\n")

    # Step 3: Calculate 3D vectors from Y to X and Y to Z.
    # We set Y as the origin (0, 0) for the map coordinates.
    # The map's y-axis (North) is opposite to the image's pixel y-axis.
    
    # Vector YX
    dx_YX = (px_X - px_Y) * m_per_pixel
    dy_YX = -(py_X - py_Y) * m_per_pixel # Negative sign to align with North
    dz_YX = h_X - h_Y
    v_YX = (dx_YX, dy_YX, dz_YX)

    # Vector YZ
    dx_YZ = (px_Z - px_Y) * m_per_pixel
    dy_YZ = -(py_Z - py_Y) * m_per_pixel # Negative sign to align with North
    dz_YZ = h_Z - h_Y
    v_YZ = (dx_YZ, dy_YZ, dz_YZ)

    print("3D Vectors on the plane (relative to point Y):")
    print(f"  Vector YX = ({v_YX[0]:.2f} m, {v_YX[1]:.2f} m, {v_YX[2]:.2f} m)")
    print(f"  Vector YZ = ({v_YZ[0]:.2f} m, {v_YZ[1]:.2f} m, {v_YZ[2]:.2f} m)\n")

    # Step 4: Find the normal vector 'n' using the cross product of v_YX and v_YZ.
    # n = v_YX x v_YZ = (ny, nx, nz)
    nx = v_YX[1] * v_YZ[2] - v_YX[2] * v_YZ[1]
    ny = v_YX[2] * v_YZ[0] - v_YX[0] * v_YZ[2]
    nz = v_YX[0] * v_YZ[1] - v_YX[1] * v_YZ[0]
    n = (nx, ny, nz)
    
    print("Normal vector 'n' to the plane (n = v_YX x v_YZ):")
    print(f"  n = ({n[0]:.2f}, {n[1]:.2f}, {n[2]:.2f})\n")

    # Step 5: Calculate the dip angle.
    # tan(dip) = sqrt(nx^2 + ny^2) / |nz|
    horizontal_component_mag = math.sqrt(n[0]**2 + n[1]**2)
    vertical_component_mag = abs(n[2])
    
    tan_dip = horizontal_component_mag / vertical_component_mag
    dip_radians = math.atan(tan_dip)
    dip_degrees = math.degrees(dip_radians)

    print("Dip Calculation:")
    print(f"  tan(dip) = sqrt(nx\u00b2 + ny\u00b2) / |nz|")
    print(f"  tan(dip) = sqrt({n[0]:.0f}\u00b2 + {n[1]:.0f}\u00b2) / |{n[2]:.0f}|")
    print(f"  tan(dip) = {horizontal_component_mag:.0f} / {vertical_component_mag:.0f}")
    print(f"  tan(dip) = {tan_dip:.4f}")
    print(f"  Dip (degrees) = atan({tan_dip:.4f}) = {dip_degrees:.2f}\u00b0\n")

    # Step 6: Round the answer to the nearest degree.
    rounded_dip = round(dip_degrees)
    print(f"Final Answer:")
    print(f"The dip of the planar surface is {dip_degrees:.2f}\u00b0, which rounds to {rounded_dip}\u00b0.")
    
    return rounded_dip

# Execute the calculation and print the final result.
final_answer = calculate_dip()
print(f"\n<<<dip = {final_answer}>>>")
