import math

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points X, Y, and Z.
    """
    # Step 1: Define the 3D coordinates of the points based on map analysis.
    # Coordinates are in the format (x_east, y_north, z_elevation) in meters.
    p_X = (545, 338, 120)
    p_Y = (195, 248, 80)
    p_Z = (445, 538, 140)

    print("Step 1: 3D Coordinates (in meters)")
    print(f"Point X: {p_X}")
    print(f"Point Y: {p_Y}")
    print(f"Point Z: {p_Z}\n")

    # Step 2: Create two vectors on the plane, using point Y as the origin.
    # Vector from Y to X
    vec_YX = (p_X[0] - p_Y[0], p_X[1] - p_Y[1], p_X[2] - p_Y[2])
    # Vector from Y to Z
    vec_YZ = (p_Z[0] - p_Y[0], p_Z[1] - p_Y[1], p_Z[2] - p_Y[2])

    print("Step 2: Vectors on the plane")
    print(f"Vector YX = X - Y = {vec_YX}")
    print(f"Vector YZ = Z - Y = {vec_YZ}\n")

    # Step 3: Calculate the normal vector N = YX x YZ.
    Nx = vec_YX[1] * vec_YZ[2] - vec_YX[2] * vec_YZ[1]
    Ny = vec_YX[2] * vec_YZ[0] - vec_YX[0] * vec_YZ[2]
    Nz = vec_YX[0] * vec_YZ[1] - vec_YX[1] * vec_YZ[0]

    print("Step 3: Normal Vector (N = YX x YZ)")
    print(f"N = ({Nx}, {Ny}, {Nz})\n")

    # Step 4: Calculate the dip angle.
    # The formula for the tangent of the dip angle is sqrt(Nx^2 + Ny^2) / |Nz|.
    tan_dip = math.sqrt(Nx**2 + Ny**2) / abs(Nz)

    print("Step 4: Calculate Dip Angle")
    print("Equation: tan(dip) = sqrt(Nx^2 + Ny^2) / |Nz|")
    print(f"tan(dip) = sqrt({Nx}^2 + {Ny}^2) / |{Nz}|")
    print(f"tan(dip) = {tan_dip:.4f}\n")

    # Calculate dip in radians and then convert to degrees.
    dip_rad = math.atan(tan_dip)
    dip_deg = math.degrees(dip_rad)

    # Round the final answer to the nearest degree.
    rounded_dip = round(dip_deg)

    print("Final Result:")
    print(f"Dip Angle = arctan({tan_dip:.4f}) = {dip_deg:.2f} degrees")
    print(f"Rounded to the nearest degree, the dip is {rounded_dip} degrees.")
    
    return rounded_dip

if __name__ == '__main__':
    dip_angle = calculate_dip()
    print(f"\n<<<9>>>")
