import numpy as np

def solve_dip_angle():
    """
    Solves the three-point problem to find the dip of a planar surface.
    """
    # Step 1: Define the data from the problem statement and map analysis.
    # Heights of the three locations
    h_X = 120.0  # meters
    h_Y = 80.0   # meters
    h_Z = 140.0  # meters

    # Pixel coordinates are measured from the image and converted to map coordinates.
    # Scale: 196 pixels = 200 meters -> scale = 200.0 / 196.0 m/pix
    # With Y at the origin, the map coordinates are:
    # X_map = (338.78, 122.45)
    # Z_map = (190.82, 317.35)
    
    # Create 3D coordinate points (x, y, z) = (Easting, Northing, Height)
    # We define the points directly to focus on the geological calculation.
    P_Y = np.array([0.0, 0.0, 80.0])
    P_X = np.array([338.78, 122.45, 120.0])
    P_Z = np.array([190.82, 317.35, 140.0])

    print("--- Step 1: Establish 3D coordinates for points X, Y, and Z ---")
    print(f"Coordinates of Y: ({P_Y[0]:.2f} m, {P_Y[1]:.2f} m, {P_Y[2]:.2f} m)")
    print(f"Coordinates of X: ({P_X[0]:.2f} m, {P_X[1]:.2f} m, {P_X[2]:.2f} m)")
    print(f"Coordinates of Z: ({P_Z[0]:.2f} m, {P_Z[1]:.2f} m, {P_Z[2]:.2f} m)\n")

    # Step 2: Create two vectors on the planar surface.
    vec_YX = P_X - P_Y
    vec_YZ = P_Z - P_Y
    
    print("--- Step 2: Define two vectors on the plane using these points ---")
    print(f"Vector from Y to X: ({vec_YX[0]:.2f}, {vec_YX[1]:.2f}, {vec_YX[2]:.2f})")
    print(f"Vector from Y to Z: ({vec_YZ[0]:.2f}, {vec_YZ[1]:.2f}, {vec_YZ[2]:.2f})\n")

    # Step 3: Calculate the normal vector to the plane.
    normal_vector = np.cross(vec_YX, vec_YZ)

    print("--- Step 3: Calculate the normal vector 'n' to the plane (via cross product) ---")
    print(f"Normal vector n = ({normal_vector[0]:.2f}, {normal_vector[1]:.2f}, {normal_vector[2]:.2f})\n")

    # Step 4: Calculate the dip angle.
    # The vertical vector is k = (0, 0, 1).
    vertical_vector = np.array([0.0, 0.0, 1.0])
    
    # Using the dot product formula to find the angle between the normal and vertical vectors.
    dot_product = np.dot(normal_vector, vertical_vector)
    norm_n = np.linalg.norm(normal_vector)
    norm_k = np.linalg.norm(vertical_vector) # This is 1
    
    # cos(angle) = |dot_product| / (norm_n * norm_k)
    cos_angle = np.abs(dot_product) / (norm_n * norm_k)
    
    dip_rad = np.arccos(cos_angle)
    dip_deg = np.degrees(dip_rad)
    
    print("--- Step 4: Calculate the dip angle ---")
    print("The dip is the angle between the normal vector 'n' and the vertical axis 'k'=(0,0,1).")
    print(f"dip = arccos( |(n . k)| / (|n|*|k|) )")
    print(f"dip = arccos( |{dot_product:.2f}| / ({norm_n:.2f} * {norm_k:.2f}) )")
    print(f"dip = arccos({cos_angle:.4f})")
    print(f"Dip in degrees = {dip_deg:.2f}°\n")

    # Step 5: Round the answer to the nearest degree.
    rounded_dip = round(dip_deg)
    
    print("--- Final Answer ---")
    print(f"The calculated dip is {dip_deg:.2f} degrees.")
    print(f"Rounding to the nearest degree, the dip of the planar surface is {rounded_dip} degrees.")
    
# Execute the function
solve_dip_angle()
print("<<<9>>>")