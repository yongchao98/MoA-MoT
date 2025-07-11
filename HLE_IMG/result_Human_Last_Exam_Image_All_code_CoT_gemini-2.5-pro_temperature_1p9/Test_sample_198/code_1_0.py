import math

def calculate_dip():
    """
    Calculates the dip of a planar surface from three points on a map.
    """
    # --- 1. Define Input Data ---
    # Heights of the locations in meters
    h_X = 120.0
    h_Y = 80.0
    h_Z = 140.0

    # Pixel coordinates measured from the map
    px_X = (1520, 690)
    px_Y = (475, 452)
    px_Z = (1031, 1400)

    # Scale measured from the map: 200 meters for 444 pixels
    map_scale = 200.0 / 444.0  # meters per pixel

    # --- 2. Convert Pixel Coordinates to Meters ---
    X_m = (px_X[0] * map_scale, px_X[1] * map_scale)
    Y_m = (px_Y[0] * map_scale, px_Y[1] * map_scale)
    Z_m = (px_Z[0] * map_scale, px_Z[1] * map_scale)

    # --- 3. Determine the Strike Line ---
    # Find a point 'W' on line YZ with the same elevation as X (120 m).
    # The line connecting X and W is the strike line.
    elevation_ratio = (h_X - h_Y) / (h_Z - h_Y)
    W_m = (
        Y_m[0] + elevation_ratio * (Z_m[0] - Y_m[0]),
        Y_m[1] + elevation_ratio * (Z_m[1] - Y_m[1])
    )

    # --- 4. Calculate Rise and Run ---
    # Rise is the vertical elevation change from the strike line to point Y.
    rise = abs(h_X - h_Y)

    # Run is the perpendicular horizontal distance from point Y to the strike line (XW).
    # Using the formula for distance from a point to a line.
    x0, y0 = Y_m
    x1, y1 = X_m
    x2, y2 = W_m
    numerator = abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1))
    denominator = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    run = numerator / denominator

    # --- 5. Calculate Dip Angle and Output Results ---
    dip_radians = math.atan(rise / run)
    dip_degrees = math.degrees(dip_radians)
    rounded_dip = round(dip_degrees)

    print("To find the dip, we use the formula: Dip = arctan(Vertical Difference / Horizontal Distance)\n")
    print("Calculations:")
    print(f"1. Vertical Difference (rise) = |Height(X) - Height(Y)| = |{h_X} m - {h_Y} m| = {rise:.2f} m")
    print(f"2. Horizontal Distance (run) perpendicular to strike = {run:.2f} m\n")

    print("Final Equation:")
    print(f"Dip = arctan({rise:.2f} / {run:.2f})")
    print(f"Dip = {dip_degrees:.1f} degrees\n")

    print(f"Rounding to the nearest degree, the dip of the planar surface is {rounded_dip} degrees.")


calculate_dip()
<<<7>>>