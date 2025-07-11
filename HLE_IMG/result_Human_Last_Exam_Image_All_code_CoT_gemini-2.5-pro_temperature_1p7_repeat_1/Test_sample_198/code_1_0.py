import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Input data extracted from the map image and problem description.
    # Elevations in meters
    h_X = 120.0
    h_Y = 80.0
    h_Z = 140.0

    # Pixel coordinates measured from the image (origin at bottom-left)
    px_X = (593, 203)
    px_Y = (235, 149)
    px_Z = (433, 442)

    # Scale bar information measured from the image
    scale_bar_start_px = 697
    scale_bar_end_px = 833
    scale_bar_dist_m = 200.0

    # Step 2: Calculate the scale factor (meters per pixel)
    scale_bar_len_px = scale_bar_end_px - scale_bar_start_px
    if scale_bar_len_px == 0:
        print("Error: Scale bar length cannot be zero.")
        return
    scale_factor = scale_bar_dist_m / scale_bar_len_px

    # Step 3: Create 3D points with real-world coordinates (x, y, z) in meters.
    P_X = (px_X[0] * scale_factor, px_X[1] * scale_factor, h_X)
    P_Y = (px_Y[0] * scale_factor, px_Y[1] * scale_factor, h_Y)
    P_Z = (px_Z[0] * scale_factor, px_Z[1] * scale_factor, h_Z)

    # Step 4: Create two vectors lying on the plane. Let's use Y as the origin.
    vec_YX = (P_X[0] - P_Y[0], P_X[1] - P_Y[1], P_X[2] - P_Y[2])
    vec_YZ = (P_Z[0] - P_Y[0], P_Z[1] - P_Y[1], P_Z[2] - P_Y[2])

    # Step 5: Calculate the normal vector to the plane using the cross product.
    # n = vec_YX x vec_YZ = (nx, ny, nz)
    nx = vec_YX[1] * vec_YZ[2] - vec_YX[2] * vec_YZ[1]
    ny = vec_YX[2] * vec_YZ[0] - vec_YX[0] * vec_YZ[2]
    nz = vec_YX[0] * vec_YZ[1] - vec_YX[1] * vec_YZ[0]

    # Step 6: Calculate the dip angle.
    # The tangent of the dip angle is the gradient of the plane.
    # tan(dip) = sqrt(nx^2 + ny^2) / |nz|
    if nz == 0:
        # This case corresponds to a vertical plane.
        dip_degrees = 90.0
        tan_dip = float('inf')
    else:
        horizontal_component = math.sqrt(nx**2 + ny**2)
        vertical_component = abs(nz)
        tan_dip = horizontal_component / vertical_component
        dip_degrees = math.degrees(math.atan(tan_dip))
    
    # Step 7: Print the calculation and the final answer.
    print("The dip angle (δ) is found using the equation: δ = arctan(gradient)")
    print("The gradient is calculated from the normal vector components (nx, ny, nz) of the plane.")
    print(f"Based on the data, the final equation is:")
    
    if nz != 0:
        print(f"Dip = arctan( {math.sqrt(nx**2 + ny**2):.1f} / |{nz:.1f}| )")
        print(f"Dip = arctan( {tan_dip:.4f} )")
        print(f"Dip = {dip_degrees:.1f}°")
    else:
        print("Dip = 90° (The plane is vertical)")
        
    rounded_dip = round(dip_degrees)
    print(f"\nRounding to the nearest degree, the dip of the planar surface is {rounded_dip} degrees.")
    
    # Final answer in specified format
    print(f"\n<<<{rounded_dip}>>>")

calculate_dip()