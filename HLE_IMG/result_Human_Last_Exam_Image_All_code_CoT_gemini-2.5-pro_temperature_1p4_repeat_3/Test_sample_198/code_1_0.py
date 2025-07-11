import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points on a map.
    """
    # Step 1: Define points and scale from map analysis.
    # Pixel coordinates (x, y) with (0,0) at top-left, measured from the image.
    X_pix = np.array([594, 409])
    Y_pix = np.array([211, 458])
    Z_pix = np.array([397, 184])

    # Heights (z-coordinate) in meters.
    h_X = 120
    h_Y = 80
    h_Z = 140

    # Scale factor (meters per pixel) from the scale bar.
    # The 100m segment of the scale bar measures 65 pixels (733 - 668).
    scale_factor = 100.0 / 65.0
    
    print("Step 1: Establishing 3D coordinates for points X, Y, and Z.")
    
    # Establish a world coordinate system with Y at the origin (0,0) on the map.
    # The x-axis points East, and the y-axis points North.
    # The image's y-pixel coordinate increases downwards (South), so we negate the change in y.
    
    # Coordinates of Y
    P_Y = np.array([0.0, 0.0, float(h_Y)])
    
    # Calculate coordinates of X relative to Y
    x_X = (X_pix[0] - Y_pix[0]) * scale_factor
    y_X = (Y_pix[1] - X_pix[1]) * scale_factor # Negating this (Y_pix[1] - X_pix[1]) makes y-axis point North
    P_X = np.array([x_X, y_X, float(h_X)])

    # Calculate coordinates of Z relative to Y
    x_Z = (Z_pix[0] - Y_pix[0]) * scale_factor
    y_Z = (Y_pix[1] - Z_pix[1]) * scale_factor
    P_Z = np.array([x_Z, y_Z, float(h_Z)])
    
    print(f"Point X: (x={P_X[0]:.1f} m, y={P_X[1]:.1f} m, z={P_X[2]:.1f} m)")
    print(f"Point Y: (x={P_Y[0]:.1f} m, y={P_Y[1]:.1f} m, z={P_Y[2]:.1f} m)")
    print(f"Point Z: (x={P_Z[0]:.1f} m, y={P_Z[1]:.1f} m, z={P_Z[2]:.1f} m)\n")

    # Step 2: Create two vectors on the plane.
    vec_YX = P_X - P_Y
    vec_YZ = P_Z - P_Y
    
    print("Step 2: Defining two vectors on the plane (from Y to X, and Y to Z).")
    print(f"Vector YX = {np.array2string(vec_YX, precision=1, separator=', ')}")
    print(f"Vector YZ = {np.array2string(vec_YZ, precision=1, separator=', ')}\n")

    # Step 3: Calculate the normal vector via cross product.
    normal_vec = np.cross(vec_YX, vec_YZ)
    nx, ny, nz = normal_vec
    
    print("Step 3: Calculating the normal vector n = (nx, ny, nz) using the cross product.")
    print(f"n = YX × YZ = {np.array2string(normal_vec, precision=1, separator=', ')}\n")

    # Step 4: Calculate the dip angle.
    print("Step 4: Calculating the dip angle (δ) using the formula: δ = arctan(sqrt(nx² + ny²) / |nz|)")
    horizontal_component_mag = np.sqrt(nx**2 + ny**2)
    tan_dip = horizontal_component_mag / np.abs(nz)
    dip_degrees = np.degrees(np.arctan(tan_dip))
    
    print(f"nx = {nx:.1f}, ny = {ny:.1f}, nz = {nz:.1f}")
    print(f"Magnitude of horizontal component = sqrt({nx:.1f}² + {ny:.1f}²) = {horizontal_component_mag:.1f}")
    print(f"Absolute value of vertical component = |{nz:.1f}| = {np.abs(nz):.1f}")
    print(f"tan(δ) = {horizontal_component_mag:.1f} / {np.abs(nz):.1f} = {tan_dip:.5f}")
    print(f"δ = arctan({tan_dip:.5f}) = {dip_degrees:.2f} degrees\n")

    # Step 5: Round to the nearest degree.
    final_dip = round(dip_degrees)
    print("Step 5: Rounding the result to the nearest degree.")
    print(f"The dip of the planar surface is {final_dip} degrees.")
    
    return final_dip

calculate_dip()