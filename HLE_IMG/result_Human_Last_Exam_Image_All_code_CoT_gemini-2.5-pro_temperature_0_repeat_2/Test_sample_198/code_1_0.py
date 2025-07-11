import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with map locations and heights.
    """
    # Step 1: Define the data from the problem description and map analysis.
    print("Step 1: Defining the known data points (locations and heights).")
    
    # Heights of the three locations in meters
    h_X = 120.0
    h_Y = 80.0
    h_Z = 140.0
    print(f"Heights: h_X = {h_X} m, h_Y = {h_Y} m, h_Z = {h_Z} m")

    # Pixel coordinates estimated from the map image.
    # (x, y) with origin at top-left corner.
    pix_X = np.array([553, 453])
    pix_Y = np.array([203, 553])
    pix_Z = np.array([403, 203])
    
    # Pixel coordinates of the scale bar to determine the scale.
    scale_0_pix = np.array([753, 753])
    scale_200_pix = np.array([903, 753])

    # Step 2: Calculate the scale (meters per pixel) and establish map coordinates.
    print("\nStep 2: Converting pixel coordinates from the map to real-world meters.")
    pixel_dist = np.linalg.norm(scale_200_pix - scale_0_pix)
    meter_dist = 200.0
    scale = meter_dist / pixel_dist
    print(f"The map scale is approximately {scale:.2f} meters per pixel.")

    # Set point Y as the origin (0, 0) of our map coordinate system.
    # The image y-axis points down, while the map y-axis (North) points up.
    # We need to invert the y-component of the relative pixel vectors.
    map_X = (pix_X - pix_Y) * scale
    map_X[1] *= -1  # Invert y-axis to align with North

    map_Z = (pix_Z - pix_Y) * scale
    map_Z[1] *= -1  # Invert y-axis to align with North

    map_Y = np.array([0.0, 0.0])
    
    print("Map coordinates relative to Y (East, North):")
    print(f"Y: ({map_Y[0]:.1f}, {map_Y[1]:.1f}) m")
    print(f"X: ({map_X[0]:.1f}, {map_X[1]:.1f}) m")
    print(f"Z: ({map_Z[0]:.1f}, {map_Z[1]:.1f}) m")

    # Step 3: Determine the equation of the planar surface, z = Ax + By + C.
    print("\nStep 3: Solving for the coefficients of the plane equation z = Ax + By + C.")
    
    # Using point Y(0, 0, 80): 80 = A*0 + B*0 + C => C = 80
    C = h_Y
    print(f"From point Y, we find C = {C:.1f}")

    # This gives a system of two linear equations for A and B:
    # 1) h_X - C = A*x_X + B*y_X
    # 2) h_Z - C = A*x_Z + B*y_Z
    M = np.array([
        [map_X[0], map_X[1]],
        [map_Z[0], map_Z[1]]
    ])
    V = np.array([h_X - C, h_Z - C])

    # Solve for A and B
    A, B = np.linalg.solve(M, V)
    print(f"Solving the system of equations gives A = {A:.4f} and B = {B:.4f}")

    # Step 4: Calculate the dip of the plane.
    print("\nStep 4: Calculating the dip angle from the plane's gradient.")
    # tan(dip) = magnitude of the gradient = sqrt(A^2 + B^2)
    tan_dip = np.sqrt(A**2 + B**2)
    
    # The dip angle in degrees is arctan(tan_dip)
    dip_degrees = np.degrees(np.arctan(tan_dip))

    print(f"The tangent of the dip angle is sqrt({A:.4f}² + {B:.4f}²) = {tan_dip:.4f}")
    print(f"The dip angle is arctan({tan_dip:.4f}) = {dip_degrees:.2f} degrees.")

    # Step 5: Round the answer to the nearest degree.
    rounded_dip = round(dip_degrees)
    print(f"\nRounding to the nearest degree, the dip is {rounded_dip} degrees.")
    
    return rounded_dip

# Run the calculation and print the final answer in the required format.
final_answer = calculate_dip()
print(f"\n<<< {final_answer} >>>")
