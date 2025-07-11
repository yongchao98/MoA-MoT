import math

def calculate_dip():
    """
    Calculates the dip of a planar surface from three points with known coordinates and heights.
    """
    # Step 1 & 2: Define coordinates and scale
    # Heights of the three locations
    h_X = 120  # meters
    h_Y = 80   # meters
    h_Z = 140  # meters

    # Pixel coordinates measured from the map image (origin at top-left)
    # Using a tool, we get: Y(195, 621), X(612, 513), Z(470, 312)
    # The scale bar from 0 to 200m is 198 pixels long (1031-833).
    scale_pixels = 198
    scale_meters = 200.0
    m_per_pixel = scale_meters / scale_pixels

    # Establish a map coordinate system with Y at the origin (0,0)
    # and Y-axis pointing North.
    # X coordinates relative to Y
    # dx is the East-West distance, dy is the North-South distance
    x_map_X = (612 - 195) * m_per_pixel
    y_map_X = (621 - 513) * m_per_pixel # Pixel y decreases going North, so we subtract X from Y.
    
    # Z coordinates relative to Y
    x_map_Z = (470 - 195) * m_per_pixel
    y_map_Z = (621 - 312) * m_per_pixel
    
    print("Step 1: Determine 3D coordinates of the points (in meters).")
    print("Point Y is set as the origin (0, 0) on the map.")
    # Step 3: Define 3D position vectors P = (x, y, h)
    P_Y = (0, 0, h_Y)
    P_X = (x_map_X, y_map_X, h_X)
    P_Z = (x_map_Z, y_map_Z, h_Z)
    print(f"Coordinates of Y: ({P_Y[0]:.2f}, {P_Y[1]:.2f}, {P_Y[2]:.2f})")
    print(f"Coordinates of X: ({P_X[0]:.2f}, {P_X[1]:.2f}, {P_X[2]:.2f})")
    print(f"Coordinates of Z: ({P_Z[0]:.2f}, {P_Z[1]:.2f}, {P_Z[2]:.2f})\n")

    # Step 4: Find the plane's normal vector using cross product
    # Vector from Y to X
    v_YX = (P_X[0] - P_Y[0], P_X[1] - P_Y[1], P_X[2] - P_Y[2])
    # Vector from Y to Z
    v_YZ = (P_Z[0] - P_Y[0], P_Z[1] - P_Y[1], P_Z[2] - P_Y[2])
    
    print("Step 2: Define two vectors on the planar surface.")
    print(f"Vector YX: ({v_YX[0]:.2f}, {v_YX[1]:.2f}, {v_YX[2]:.2f})")
    print(f"Vector YZ: ({v_YZ[0]:.2f}, {v_YZ[1]:.2f}, {v_YZ[2]:.2f})\n")

    # Normal vector N = v_YX x v_YZ = (nx, ny, nh)
    nx = v_YX[1] * v_YZ[2] - v_YX[2] * v_YZ[1]
    ny = v_YX[2] * v_YZ[0] - v_YX[0] * v_YZ[2]
    nh = v_YX[0] * v_YZ[1] - v_YX[1] * v_YZ[0]
    
    print("Step 3: Calculate the normal vector N to the plane using the cross product of YX and YZ.")
    print(f"Normal vector N = (nx, ny, nh) = ({nx:.2f}, {ny:.2f}, {nh:.2f})\n")

    # Step 5: Calculate the dip angle
    # tan(dip) = sqrt(nx^2 + ny^2) / |nh|
    tan_dip = math.sqrt(nx**2 + ny**2) / abs(nh)
    dip_radians = math.atan(tan_dip)
    dip_degrees = math.degrees(dip_radians)

    print("Step 4: Calculate the dip using the formula tan(dip) = sqrt(nx^2 + ny^2) / |nh|")
    print(f"tan(dip) = sqrt({nx:.2f}^2 + {ny:.2f}^2) / |{nh:.2f}|")
    print(f"tan(dip) = {tan_dip:.4f}")
    print(f"Dip = atan({tan_dip:.4f}) = {dip_degrees:.2f} degrees\n")

    # Step 6: Final Result
    dip_rounded = round(dip_degrees)
    print("Step 5: Round the result to the nearest degree.")
    print(f"The dip of the planar surface is {dip_rounded} degrees.")
    
    return dip_rounded

# Run the calculation and store the final answer
final_answer = calculate_dip()
# The final answer is also printed inside the function.
# To match the requested final format, we can capture it this way.
print(f"<<<{final_answer}>>>")