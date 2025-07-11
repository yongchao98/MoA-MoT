import math

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points X, Y, and Z.
    """
    # --- Step 1: Establish Coordinate System ---
    # Heights from the problem description
    h_X = 120  # meters
    h_Y = 80   # meters
    h_Z = 140  # meters

    # Pixel coordinates measured from the map image (e.g., using an image editor)
    # Origin (0,0) of pixel measurement is top-left corner of the image.
    px_X = (585, 517)
    px_Y = (195, 608)
    px_Z = (447, 280)

    # Scale bar measurement: 200 meters corresponds to 208 pixels in the image.
    # (Pixel at 200m mark - Pixel at 0m mark = 1008 - 800 = 208 pixels)
    scale_m_per_px = 200.0 / 208.0

    # Calculate map coordinates in meters, with Y as the origin (0, 0).
    # The map's y-axis points North (up), while pixel y-coordinates increase downwards.
    # So a smaller pixel y-coordinate means a larger North coordinate.
    x_Y, y_Y = 0.0, 0.0
    x_X = (px_X[0] - px_Y[0]) * scale_m_per_px
    y_X = (px_Y[1] - px_X[1]) * scale_m_per_px # (px_Y > px_X -> Y is south of X)
    x_Z = (px_Z[0] - px_Y[0]) * scale_m_per_px
    y_Z = (px_Y[1] - px_Z[1]) * scale_m_per_px # (px_Y > px_Z -> Y is south of Z)
    
    print("--- Calculated Map and 3D Coordinates (in meters) ---")
    print(f"Point X: (x={x_X:.1f}, y={y_X:.1f}, z={h_X})")
    print(f"Point Y: (x={x_Y:.1f}, y={y_Y:.1f}, z={h_Y})")
    print(f"Point Z: (x={x_Z:.1f}, y={y_Z:.1f}, z={h_Z})")
    print("-" * 50)

    # --- Step 2 & 3: Define 3D points and vectors on the plane ---
    # Vector from Y to X
    vec_YX = (x_X - x_Y, y_X - y_Y, h_X - h_Y)
    # Vector from Y to Z
    vec_YZ = (x_Z - x_Y, y_Z - y_Y, h_Z - h_Y)
    
    print("--- Vectors on the Plane ---")
    print(f"Vector YX = ({vec_YX[0]:.1f}, {vec_YX[1]:.1f}, {vec_YX[2]:.1f})")
    print(f"Vector YZ = ({vec_YZ[0]:.1f}, {vec_YZ[1]:.1f}, {vec_YZ[2]:.1f})")
    print("-" * 50)
    
    # --- Step 4: Find the Normal Vector (Cross Product) ---
    # n = YX x YZ
    nx = vec_YX[1] * vec_YZ[2] - vec_YX[2] * vec_YZ[1]
    ny = vec_YX[2] * vec_YZ[0] - vec_YX[0] * vec_YZ[2]
    nz = vec_YX[0] * vec_YZ[1] - vec_YX[1] * vec_YZ[0]
    
    normal_vector = (nx, ny, nz)
    
    print("--- Normal Vector Calculation ---")
    print(f"Normal vector n = (nx, ny, nz) = ({nx:.1f}, {ny:.1f}, {nz:.1f})")
    print("-" * 50)

    # --- Step 5: Calculate the Dip Angle ---
    # Magnitude of the normal vector
    norm_n = math.sqrt(nx**2 + ny**2 + nz**2)
    
    # The z-component of the normal vector
    abs_nz = abs(nz)

    # Cosine of the dip angle
    cos_delta = abs_nz / norm_n

    # Dip angle in radians
    delta_rad = math.acos(cos_delta)
    
    # Dip angle in degrees
    delta_deg = math.degrees(delta_rad)
    
    # Rounded dip angle
    rounded_dip = round(delta_deg)
    
    print("--- Dip Angle Calculation ---")
    print("The dip angle (δ) is calculated using the formula: δ = arccos(|nz| / ||n||)")
    print(f"The equation with values is: δ = arccos(|{nz:.1f}| / {norm_n:.1f})")
    print(f"|{nz:.1f}| / {norm_n:.1f} = {cos_delta:.4f}")
    print(f"δ = arccos({cos_delta:.4f}) = {delta_rad:.4f} radians")
    print(f"In degrees, the dip angle is {delta_deg:.2f}°")
    print("-" * 50)
    
    print(f"The dip of the planar surface, rounded to the nearest degree, is {rounded_dip} degrees.")
    
    return rounded_dip

if __name__ == '__main__':
    final_answer = calculate_dip()
    print(f"<<<{final_answer}>>>")