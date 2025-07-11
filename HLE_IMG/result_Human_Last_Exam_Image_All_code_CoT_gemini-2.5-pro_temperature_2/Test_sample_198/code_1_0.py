import numpy as np
import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with heights on a map.
    """
    # Step 1: Define the data from the problem
    # Heights of the locations
    h_X = 120  # meters
    h_Y = 80   # meters
    h_Z = 140  # meters

    # Pixel coordinates measured from the image (origin at top-left)
    # These values are determined by inspecting the image file.
    px_X = (500, 396)
    px_Y = (252, 461)
    px_Z = (350, 201)

    # Pixel coordinates for the scale bar
    scale_bar_start_px = 590
    scale_bar_end_px = 720
    scale_bar_dist_m = 200

    print("Step 1: Analyzing the map data")
    print(f"Height at X: {h_X} m, Height at Y: {h_Y} m, Height at Z: {h_Z} m")
    print("-" * 30)

    # Step 2: Calculate the map scale
    pixel_dist_scale = scale_bar_end_px - scale_bar_start_px
    scale_m_per_px = scale_bar_dist_m / pixel_dist_scale

    print("Step 2: Determining the map scale")
    print(f"The 200m scale bar measures {pixel_dist_scale} pixels.")
    print(f"Therefore, the scale is {scale_m_per_px:.4f} meters per pixel.")
    print("-" * 30)

    # Step 3: Establish a map coordinate system and find coordinates of X and Z relative to Y
    # Set Y as the origin (0, 0). East is +x, North is +y.
    # Pixel y-axis increases downwards, so a negative sign is needed for Northing.
    rel_X_map_x = (px_X[0] - px_Y[0]) * scale_m_per_px
    rel_X_map_y = -(px_X[1] - px_Y[1]) * scale_m_per_px

    rel_Z_map_x = (px_Z[0] - px_Y[0]) * scale_m_per_px
    rel_Z_map_y = -(px_Z[1] - px_Y[1]) * scale_m_per_px

    # Height differences relative to Y
    delta_h_X = h_X - h_Y
    delta_h_Z = h_Z - h_Y
    
    print("Step 3: Calculating relative positions")
    print("Using point Y as the origin (0,0) with North as the positive y-axis:")
    print(f"Position of X relative to Y: ({rel_X_map_x:.1f} m East, {rel_X_map_y:.1f} m North)")
    print(f"Position of Z relative to Y: ({rel_Z_map_x:.1f} m East, {rel_Z_map_y:.1f} m North)")
    print(f"Height of X relative to Y: {delta_h_X} m")
    print(f"Height of Z relative to Y: {delta_h_Z} m")
    print("-" * 30)

    # Step 4: Solve for the plane's gradient components a and b
    # The plane is defined by Δh = a*x + b*y. We have a system of two linear equations:
    # A * [a, b]' = B
    A = np.array([
        [rel_X_map_x, rel_X_map_y],
        [rel_Z_map_x, rel_Z_map_y]
    ])
    B = np.array([delta_h_X, delta_h_Z])

    try:
        gradient_vector = np.linalg.solve(A, B)
        a, b = gradient_vector
    except np.linalg.LinAlgError:
        print("Error: The points are collinear, and a unique plane cannot be determined.")
        return

    print("Step 4: Solving for the plane's gradient (a, b)")
    print("From the linear equations:")
    print(f"  {delta_h_X} = a * {rel_X_map_x:.1f} + b * {rel_X_map_y:.1f}")
    print(f"  {delta_h_Z} = a * {rel_Z_map_x:.1f} + b * {rel_Z_map_y:.1f}")
    print(f"The solved gradient components are: a = {a:.4f}, b = {b:.4f}")
    print("-" * 30)

    # Step 5: Calculate the dip
    # The tangent of the dip angle is the magnitude of the gradient vector.
    tan_dip = math.sqrt(a**2 + b**2)
    dip_rad = math.atan(tan_dip)
    dip_deg = math.degrees(dip_rad)
    rounded_dip = round(dip_deg)

    print("Step 5: Calculating the dip from the gradient")
    print("The final equation for the dip angle (δ) is: δ = arctan(sqrt(a² + b²))")
    print("\nPlugging in the calculated values:")
    print(f"tan(δ) = sqrt({a:.4f}² + {b:.4f}²)")
    print(f"tan(δ) = sqrt({a**2:.5f} + {b**2:.5f})")
    print(f"tan(δ) = {tan_dip:.4f}")
    print(f"δ = arctan({tan_dip:.4f}) = {dip_deg:.2f}°")
    print("\n" + "=" * 30)
    print(f"The calculated dip is {dip_deg:.2f} degrees.")
    print(f"Rounded to the nearest degree, the dip is {rounded_dip} degrees.")
    print("=" * 30)
    return rounded_dip

if __name__ == '__main__':
    final_answer = calculate_dip()
    print(f"\n<<<>>>")
    # The line above with <<<>>> is a placeholder. 
    # In a real application, you might use it to return the value.
    # For this script, we'll print it explicitly to be clear.
    print(f"<<<{final_answer}>>>")