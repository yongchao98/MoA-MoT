import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with their
    heights and map positions.
    """
    # Step 1: Define the given data.
    # Heights of the three points in meters.
    h_X = 120
    h_Y = 80
    h_Z = 140

    print("--- Input Data ---")
    print(f"Height of point X: {h_X} m")
    print(f"Height of point Y: {h_Y} m")
    print(f"Height of point Z: {h_Z} m")
    print("\n")

    # Step 2: Determine map coordinates from the image.
    # These values are derived from analyzing the map image and its scale.
    # The map's scale is 200 meters for 188 pixels.
    scale = 200.0 / 188.0  # meters per pixel

    # We set point Y as the origin (0, 0) of our map's (x, y) coordinate system.
    # The North arrow is up, so the map's y-axis points North.
    # Pixel coordinates relative to Y are measured from the image:
    # X relative to Y: (482, -58) pixels (dx, dy). dy is negative as image y-axis is inverted.
    # Z relative to Y: (192, -306) pixels (dx, dy).
    map_x_X = 482 * scale
    map_y_X = -(-58) * scale # Convert to North
    map_x_Z = 192 * scale
    map_y_Z = -(-306) * scale # Convert to North

    # Combine map coordinates with heights to get 3D points.
    p_X = (map_x_X, map_y_X, float(h_X))
    p_Y = (0.0, 0.0, float(h_Y))
    p_Z = (map_x_Z, map_y_Z, float(h_Z))

    print("--- 3D Coordinates (Y at origin) ---")
    print(f"Point X: ({p_X[0]:.2f}, {p_X[1]:.2f}, {p_X[2]:.2f}) m")
    print(f"Point Y: ({p_Y[0]:.2f}, {p_Y[1]:.2f}, {p_Y[2]:.2f}) m")
    print(f"Point Z: ({p_Z[0]:.2f}, {p_Z[1]:.2f}, {p_Z[2]:.2f}) m")
    print("\n")

    # Step 3: Define two vectors on the plane.
    vec_YX = (p_X[0] - p_Y[0], p_X[1] - p_Y[1], p_X[2] - p_Y[2])
    vec_YZ = (p_Z[0] - p_Y[0], p_Z[1] - p_Y[1], p_Z[2] - p_Y[2])
    
    # Step 4: Calculate the normal vector (N) using the cross product.
    # N = vec_YX x vec_YZ
    N_x = vec_YX[1] * vec_YZ[2] - vec_YX[2] * vec_YZ[1]
    N_y = vec_YX[2] * vec_YZ[0] - vec_YX[0] * vec_YZ[2]
    N_z = vec_YX[0] * vec_YZ[1] - vec_YX[1] * vec_YZ[0]
    
    print("--- Normal Vector Calculation ---")
    print(f"Vector YX = <{vec_YX[0]:.2f}, {vec_YX[1]:.2f}, {vec_YX[2]:.2f}>")
    print(f"Vector YZ = <{vec_YZ[0]:.2f}, {vec_YZ[1]:.2f}, {vec_YZ[2]:.2f}>")
    print(f"Normal Vector N = YX x YZ = <{N_x:.2f}, {N_y:.2f}, {N_z:.2f}>")
    print("\n")

    # Step 5: Calculate the dip angle (δ).
    # tan(δ) = |horizontal component of N| / |vertical component of N|
    horizontal_magnitude = math.sqrt(N_x**2 + N_y**2)
    vertical_magnitude = abs(N_z)
    
    tan_dip = horizontal_magnitude / vertical_magnitude
    dip_radians = math.atan(tan_dip)
    dip_degrees = math.degrees(dip_radians)

    # Step 6: Round the answer to the nearest degree.
    rounded_dip = round(dip_degrees)

    print("--- Dip Calculation ---")
    print(f"The equation for the tangent of the dip angle (δ) is: tan(δ) = sqrt(Nx² + Ny²) / |Nz|")
    print(f"Plugging in the numbers:")
    print(f"tan(δ) = sqrt({N_x:.2f}² + {N_y:.2f}²) / |{N_z:.2f}|")
    print(f"tan(δ) = {horizontal_magnitude:.2f} / {vertical_magnitude:.2f} = {tan_dip:.5f}")
    print(f"δ = arctan({tan_dip:.5f}) = {dip_degrees:.2f} degrees")
    print("\n")

    print("--- Final Answer ---")
    print(f"The dip of the planar surface rounded to the nearest degree is: {rounded_dip} degrees.")
    
    return rounded_dip

# Execute the calculation and print the final answer in the required format.
final_answer = calculate_dip()
# print(f"\n<<<{final_answer}>>>") # This is commented out to avoid duplicate final answer display
# The final answer will be returned by the function and wrapped by the calling environment.

calculate_dip()
<<<9>>>