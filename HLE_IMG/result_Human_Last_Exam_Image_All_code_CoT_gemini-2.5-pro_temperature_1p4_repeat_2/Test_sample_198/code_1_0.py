import math
import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with elevations on a map.
    """
    # Step 1: Input data from map analysis.
    # Coordinates are measured in pixels from the top-left corner of the image.
    # Heights are in meters.
    X_pix = np.array([460, 375])
    h_X = 120.0
    Y_pix = np.array([138, 423])
    h_Y = 80.0
    Z_pix = np.array([320, 185])
    h_Z = 140.0

    # Scale bar information in pixels and meters.
    scale_bar_start_pix = 613
    scale_bar_end_pix = 753
    scale_bar_dist_m = 200.0

    # Step 2: Calculate the scale factor (meters per pixel).
    scale_factor = scale_bar_dist_m / (scale_bar_end_pix - scale_bar_start_pix)

    # Step 3: Determine a strike line.
    # Find a point P on the line segment YZ that has the same elevation as X (120 m).
    # We use linear interpolation based on height.
    # The fractional distance of P along the segment YZ is given by the ratio of height differences.
    ratio = (h_X - h_Y) / (h_Z - h_Y)
    
    # Calculate the pixel coordinates of point P.
    # P = Y + ratio * (Z - Y)
    vec_YZ = Z_pix - Y_pix
    P_pix = Y_pix + ratio * vec_YZ

    # Step 4: Calculate the horizontal 'run' and vertical 'rise'.
    # The 'rise' is the elevation difference between the strike line (h_X) and point Y (h_Y).
    vertical_rise_m = abs(h_X - h_Y)

    # The 'run' is the perpendicular distance from point Y to the strike line (XP).
    # We can calculate this using the vector cross product formula for distance from a point to a line:
    # Distance = |cross(vec_XP, vec_YP)| / |magnitude(vec_XP)|
    vec_XP = X_pix - P_pix
    vec_YP = Y_pix - P_pix
    
    # For 2D vectors, the magnitude of the cross product is |x1*y2 - y1*x2|.
    cross_product_mag = np.abs(np.cross(vec_XP, vec_YP))
    magnitude_XP = np.linalg.norm(vec_XP)
    
    horizontal_dist_pix = cross_product_mag / magnitude_XP
    horizontal_run_m = horizontal_dist_pix * scale_factor

    # Step 5: Calculate the dip.
    # The dip angle's tangent is the ratio of rise to run.
    tan_dip = vertical_rise_m / horizontal_run_m
    dip_rad = math.atan(tan_dip)
    dip_deg = math.degrees(dip_rad)
    rounded_dip = int(round(dip_deg))

    # Output the final results, including the equation with values.
    print("To find the dip, we use the formula: tan(dip) = Vertical Rise / Horizontal Run")
    print(f"The Vertical Rise is the elevation difference: |{h_X:.0f} m - {h_Y:.0f} m| = {vertical_rise_m:.0f} m.")
    print(f"The Horizontal Run is the map distance from point Y to the 120m strike line: {horizontal_run_m:.2f} m.")
    print("\nPlugging these numbers into the equation:")
    print(f"tan(dip) = {vertical_rise_m:.0f} / {horizontal_run_m:.2f}")
    print(f"dip = arctan({tan_dip:.4f}) = {dip_deg:.2f} degrees")
    print("\n-------------------------------------------------")
    print(f"The final answer, rounded to the nearest degree, is {rounded_dip} degrees.")
    print("-------------------------------------------------")


calculate_dip()
<<<8>>>