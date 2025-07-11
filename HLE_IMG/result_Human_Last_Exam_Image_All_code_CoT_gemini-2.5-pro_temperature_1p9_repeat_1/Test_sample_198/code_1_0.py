import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with heights.
    """

    # Step 1: Define the 3D coordinates of the points.
    # We place point Y at the origin (0,0) of the map.
    # By observing the map and using the scale bar, we estimate the coordinates for X and Z.
    # The map coordinates are (East, North) and heights are elevation.
    # Points are in the format (x, y, z) in meters.
    
    H_X = 120.0  # Height at X
    H_Y = 80.0   # Height at Y
    H_Z = 140.0  # Height at Z
    
    # Based on the map, we estimate the coordinates of X and Z relative to Y.
    # Let Y = (0, 0) on the map.
    # X appears to be about 240 m East and 80 m North of Y.
    # Z appears to be about 120 m East and 240 m North of Y.
    P_Y = np.array([0.0, 0.0, H_Y])
    P_X = np.array([240.0, 80.0, H_X])
    P_Z = np.array([120.0, 240.0, H_Z])

    print("Step 1: Define the 3D coordinates of the points (in meters).")
    print(f"Point Y: ({P_Y[0]}, {P_Y[1]}, {P_Y[2]})")
    print(f"Point X: ({P_X[0]}, {P_X[1]}, {P_X[2]})")
    print(f"Point Z: ({P_Z[0]}, {P_Z[1]}, {P_Z[2]})")
    print("-" * 30)

    # Step 2: Calculate two vectors on the plane.
    v_YX = P_X - P_Y
    v_YZ = P_Z - P_Y

    print("Step 2: Define two vectors on the planar surface from point Y.")
    print(f"Vector YX = X - Y = ({v_YX[0]}, {v_YX[1]}, {v_YX[2]})")
    print(f"Vector YZ = Z - Y = ({v_YZ[0]}, {v_YZ[1]}, {v_YZ[2]})")
    print("-" * 30)

    # Step 3: Calculate the normal vector N = (A, B, C) using the cross product.
    N = np.cross(v_YX, v_YZ)
    A, B, C = N[0], N[1], N[2]
    
    print("Step 3: Find the normal vector N = (A, B, C) to the plane.")
    print("N = YX x YZ")
    print(f"N = ({A}, {B}, {C})")
    print("-" * 30)

    # Step 4: Calculate the dip angle (delta)
    print("Step 4: Calculate the dip angle (δ) using the normal vector components.")
    print("The formula is: tan(δ) = sqrt(A² + B²) / |C|")
    
    A2_plus_B2 = A**2 + B**2
    sqrt_A2_plus_B2 = np.sqrt(A2_plus_B2)
    abs_C = abs(C)
    tan_dip = sqrt_A2_plus_B2 / abs_C
    dip_rad = np.arctan(tan_dip)
    dip_deg = np.degrees(dip_rad)
    
    print(f"tan(δ) = sqrt(({A})² + ({B})²) / |{C}|")
    print(f"tan(δ) = sqrt({A**2} + {B**2}) / {abs_C}")
    print(f"tan(δ) = sqrt({A2_plus_B2}) / {abs_C}")
    print(f"tan(δ) = {sqrt_A2_plus_B2:.4f} / {abs_C:.4f}")
    print(f"tan(δ) = {tan_dip:.4f}")
    print(f"δ = arctan({tan_dip:.4f}) = {dip_deg:.2f}°")
    print("-" * 30)
    
    # Step 5: Round to the nearest degree
    dip_rounded = round(dip_deg)
    print("Step 5: Final Answer")
    print(f"The calculated dip is {dip_deg:.2f} degrees.")
    print(f"Rounded to the nearest degree, the dip of the planar surface is: {dip_rounded}°")
    
    return dip_rounded

if __name__ == '__main__':
    final_answer = calculate_dip()
    # The final answer is enclosed in <<<>>>
    print(f"\n<<< {final_answer} >>>")