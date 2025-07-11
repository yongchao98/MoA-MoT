import math

def calculate_dip():
    """
    Calculates the dip of a planar surface from three points on a map.
    """
    # Step 1: Define initial data from the problem description and map analysis
    # Heights in meters
    heights = {'X': 120, 'Y': 80, 'Z': 140}

    # Pixel coordinates estimated from the map (origin at top-left)
    # The exact pixel values might vary slightly, but relative positions are key.
    # The image is 1000x800 px. We will flip the y-axis for standard Cartesian coordinates.
    image_height = 800
    coords_px = {
        'X': (515, image_height - 435),
        'Y': (170, image_height - 545),
        'Z': (390, image_height - 240)
    }
    
    # Scale bar pixel coordinates
    scale_bar_start_px = (725, image_height - 715)
    scale_bar_end_px = (895, image_height - 715)
    scale_bar_meters = 200

    # Step 2: Calculate the scale (meters per pixel)
    scale_bar_pixels = math.sqrt((scale_bar_end_px[0] - scale_bar_start_px[0])**2 + 
                                 (scale_bar_end_px[1] - scale_bar_start_px[1])**2)
    scale = scale_bar_meters / scale_bar_pixels

    # Convert pixel coordinates to map coordinates in meters
    coords_m = {p: (c[0] * scale, c[1] * scale) for p, c in coords_px.items()}

    # Step 3: Find point W on line YZ with the same elevation as X (120m)
    # The elevation change from Y to Z is linear along the line segment YZ.
    h_X, h_Y, h_Z = heights['X'], heights['Y'], heights['Z']
    
    # The ratio of the distance YW to YZ is equal to the ratio of elevation change
    ratio = (h_X - h_Y) / (h_Z - h_Y)
    
    # Calculate coordinates of W using vector addition: W = Y + ratio * (Z - Y)
    Y_m, Z_m = coords_m['Y'], coords_m['Z']
    W_m = (
        Y_m[0] + ratio * (Z_m[0] - Y_m[0]),
        Y_m[1] + ratio * (Z_m[1] - Y_m[1])
    )
    
    # Line XW is the strike line at 120m elevation.

    # Step 4: Calculate the 'run' and 'rise'
    # The 'run' is the perpendicular horizontal distance from point Y to the strike line XW.
    # We can calculate this using the area of the parallelogram formed by vectors XY and XW.
    # run = |Area| / |base| = |(XW) x (XY)| / |XW|
    X_m = coords_m['X']
    
    # Vector XW
    v_XW = (W_m[0] - X_m[0], W_m[1] - X_m[1])
    # Vector XY
    v_XY = (Y_m[0] - X_m[0], Y_m[1] - X_m[1])
    
    # Magnitude of 2D cross product |v1.x*v2.y - v1.y*v2.x|
    cross_product_mag = abs(v_XW[0] * v_XY[1] - v_XW[1] * v_XY[0])
    
    # Magnitude of vector XW
    mag_v_XW = math.sqrt(v_XW[0]**2 + v_XW[1]**2)
    
    run = cross_product_mag / mag_v_XW
    
    # The 'rise' is the vertical elevation difference between point Y and the strike line.
    rise = h_X - h_Y # Elevation of strike line (120m) - Elevation of Y (80m)

    # Step 5: Calculate the dip angle
    dip_radians = math.atan(rise / run)
    dip_degrees = math.degrees(dip_radians)
    rounded_dip = round(dip_degrees)

    print("Step-by-step calculation of the dip angle:")
    print(f"1. A strike line is defined by connecting point X (elevation {heights['X']}m) with a calculated point W on the line YZ which also has an elevation of {heights['X']}m.")
    print(f"2. The 'rise' is the elevation difference between the strike line and point Y: {heights['X']}m - {heights['Y']}m = {rise:.2f} m.")
    print(f"3. The 'run' is the horizontal distance from point Y to the strike line: {run:.2f} m.")
    print(f"4. The dip angle is calculated using the formula: dip = arctan(rise / run).")
    print(f"   dip = arctan({rise:.2f} / {run:.2f})")
    print(f"   dip = arctan({rise / run:.4f})")
    print(f"   dip = {dip_degrees:.2f} degrees")
    print(f"\nRounding to the nearest degree, the dip of the planar surface is: {rounded_dip} degrees.")
    
    return rounded_dip

final_answer = calculate_dip()
print(f"\n<<< {final_answer} >>>")
