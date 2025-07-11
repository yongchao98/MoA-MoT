import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points X, Y, and Z.
    """
    # Step 1 & 2: Define 3D points based on map and height data.
    
    # Heights from the problem description
    h_x = 120  # meters
    h_y = 80   # meters
    h_z = 140  # meters

    # Pixel coordinates measured from the image (origin at top-left).
    # North is up (decreasing y_pixel).
    # Y_pix = (x, y)
    Y_pix = (252, 574)
    X_pix = (690, 484)
    Z_pix = (416, 260)
    
    # Scale from the map: 200 meters corresponds to 196 pixels.
    scale = 200.0 / 196.0  # meters per pixel

    # Set point Y as the origin of our 3D coordinate system.
    Y = np.array([0.0, 0.0, h_y])

    # Calculate map coordinates for X and Z relative to Y.
    # Easting (x) is horizontal on the map. Northing (y) is vertical.
    x_coord_X = (X_pix[0] - Y_pix[0]) * scale
    y_coord_X = -(X_pix[1] - Y_pix[1]) * scale # Negative because y_pixel decreases going North
    X = np.array([x_coord_X, y_coord_X, h_x])

    x_coord_Z = (Z_pix[0] - Y_pix[0]) * scale
    y_coord_Z = -(Z_pix[1] - Y_pix[1]) * scale # Negative because y_pixel decreases going North
    Z = np.array([x_coord_Z, y_coord_Z, h_z])

    print("Step 1 & 2: 3D Coordinates (in meters, Y is origin)")
    print(f"Point Y: ({Y[0]:.2f}, {Y[1]:.2f}, {Y[2]:.2f})")
    print(f"Point X: ({X[0]:.2f}, {X[1]:.2f}, {X[2]:.2f})")
    print(f"Point Z: ({Z[0]:.2f}, {Z[1]:.2f}, {Z[2]:.2f})\n")

    # Step 3: Find the plane's normal vector.
    # Create two vectors on the plane.
    vec_YX = X - Y
    vec_YZ = Z - Y
    
    # Calculate the cross product to find the normal vector N = <a, b, c>.
    N = np.cross(vec_YX, vec_YZ)
    a, b, c = N[0], N[1], N[2]
    
    print("Step 3: Calculate the Normal Vector")
    print(f"Vector YX = <{vec_YX[0]:.2f}, {vec_YX[1]:.2f}, {vec_YX[2]:.2f}>")
    print(f"Vector YZ = <{vec_YZ[0]:.2f}, {vec_YZ[1]:.2f}, {vec_YZ[2]:.2f}>")
    print(f"Normal Vector N = <a, b, c> = <{a:.2f}, {b:.2f}, {c:.2f}>\n")

    # Step 4: Calculate the dip angle.
    print("Step 4: Calculate the Dip Angle")
    # tan(dip) = sqrt(a^2 + b^2) / |c|
    sqrt_a2_b2 = np.sqrt(a**2 + b**2)
    abs_c = np.abs(c)
    tan_dip = sqrt_a2_b2 / abs_c
    
    dip_rad = np.arctan(tan_dip)
    dip_deg = np.degrees(dip_rad)
    
    # Print the equation with the calculated values
    print("The formula for the tangent of the dip angle is: tan(dip) = sqrt(a² + b²) / |c|")
    print("Plugging in the numbers from the normal vector:")
    print(f"tan(dip) = sqrt(({a:.2f})² + ({b:.2f})²) / |{c:.2f}|")
    print(f"tan(dip) = {sqrt_a2_b2:.2f} / {abs_c:.2f}")
    print(f"tan(dip) = {tan_dip:.4f}\n")
    
    # Step 5: Final Result
    print("Step 5: Final Result")
    print(f"Dip = arctan({tan_dip:.4f})")
    print(f"Dip = {dip_deg:.2f} degrees")
    
    rounded_dip = round(dip_deg)
    print(f"\nThe dip of the planar surface rounded to the nearest degree is: {rounded_dip} degrees.")
    return rounded_dip

# Execute the function and get the final answer.
final_answer = calculate_dip()
print(f"\n<<< {final_answer} >>>")
