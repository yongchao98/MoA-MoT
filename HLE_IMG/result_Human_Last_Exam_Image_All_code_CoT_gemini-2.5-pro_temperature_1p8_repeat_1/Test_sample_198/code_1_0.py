import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points on a map.
    """
    # Step 1: Define points and their properties based on the map.
    # Heights of the three locations in meters.
    h_X = 120
    h_Y = 80
    h_Z = 140

    # Pixel coordinates are measured from the image with the origin at the top-left.
    # An image analysis tool was used to obtain these values.
    X_px = (558, 388)
    Y_px = (211, 258)
    Z_px = (454, 608)

    # The scale bar shows that 200 meters correspond to 198 pixels in the image.
    scale_px = 198.0
    scale_m = 200.0
    # Calculate the conversion factor from pixels to meters.
    pixels_to_meters_ratio = scale_m / scale_px

    print("Step 1: Determine the 3D coordinates of points X, Y, and Z.")
    # We establish a 3D coordinate system (x=East, y=North, z=Height).
    # The map's North direction corresponds to the negative pixel y-axis.
    # We'll place point Y at the map origin (0,0) for simplicity.

    # Calculate coordinates of X and Z in meters, relative to Y.
    # x is East (positive pixel x), y is North (negative pixel y).
    x_X = (X_px[0] - Y_px[0]) * pixels_to_meters_ratio
    y_X = -(X_px[1] - Y_px[1]) * pixels_to_meters_ratio
    x_Z = (Z_px[0] - Y_px[0]) * pixels_to_meters_ratio
    y_Z = -(Z_px[1] - Y_px[1]) * pixels_to_meters_ratio

    # Define the 3D points in space.
    P_Y = np.array([0.0, 0.0, float(h_Y)])
    P_X = np.array([x_X, y_X, float(h_X)])
    P_Z = np.array([x_Z, y_Z, float(h_Z)])
    
    print(f"The scale is {scale_m} meters per {scale_px} pixels.")
    print(f"With Y as the origin (0, 0, {h_Y}), the 3D coordinates are:")
    print(f"X = ({P_X[0]:.1f} m, {P_X[1]:.1f} m, {P_X[2]:.1f} m)")
    print(f"Z = ({P_Z[0]:.1f} m, {P_Z[1]:.1f} m, {P_Z[2]:.1f} m)\n")

    # Step 2: Create two vectors lying on the plane.
    v_YX = P_X - P_Y
    v_YZ = P_Z - P_Y

    print("Step 2: Create two vectors on the planar surface from these points.")
    print(f"Vector YX = <{v_YX[0]:.1f}, {v_YX[1]:.1f}, {v_YX[2]:.1f}>")
    print(f"Vector YZ = <{v_YZ[0]:.1f}, {v_YZ[1]:.1f}, {v_YZ[2]:.1f}>\n")

    # Step 3: Find the normal vector to the plane via the cross product.
    normal_vector = np.cross(v_YX, v_YZ)
    nx, ny, nz = normal_vector

    print("Step 3: Calculate the normal vector (n) to the plane using the cross product of YX and YZ.")
    print(f"n = <nx, ny, nz> = <{nx:.1f}, {ny:.1f}, {nz:.1f}>\n")

    # Step 4: Calculate the dip angle.
    # The tangent of the dip angle is the magnitude of the plane's gradient.
    # tan(dip) = sqrt(nx^2 + ny^2) / |nz|
    tan_dip = np.sqrt(nx**2 + ny**2) / np.abs(nz)
    dip_radians = np.arctan(tan_dip)
    dip_degrees = np.degrees(dip_radians)
    dip_rounded = round(dip_degrees)
    
    print("Step 4: Calculate the dip angle (δ) using the components of the normal vector.")
    print("The formula is: δ = arctan(sqrt(nx² + ny²) / |nz|)")
    print(f"δ = arctan(sqrt({nx:.1f}² + ({ny:.1f})²) / |{nz:.1f}|)")
    print(f"δ = arctan(sqrt({nx**2:.1f} + {ny**2:.1f}) / {np.abs(nz):.1f})")
    print(f"δ = arctan({np.sqrt(nx**2 + ny**2):.1f} / {np.abs(nz):.1f})")
    print(f"δ = arctan({tan_dip:.4f})")
    print(f"δ = {dip_degrees:.2f} degrees\n")

    print(f"The calculated dip of the planar surface is {dip_degrees:.2f} degrees.")
    print(f"Rounding to the nearest degree, the dip is {dip_rounded} degrees.")
    
    # Return the final rounded answer.
    return dip_rounded

# Run the calculation and print the final answer in the specified format.
final_answer = calculate_dip()
print(f"<<<{final_answer}>>>")
