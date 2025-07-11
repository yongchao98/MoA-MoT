import math

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define input data from the map and problem description.
    # Pixel coordinates (x, y) were measured from the image (origin at top-left).
    px_X = (612, 455)
    px_Y = (211, 303)
    px_Z = (438, 626)

    # Heights (elevations) in meters.
    h_X = 120
    h_Y = 80
    h_Z = 140

    # Scale bar information (in pixels).
    scale_0_px = 787
    scale_200_px = 985
    scale_dist_m = 200

    print("Step 1: Analyzing the map data...")
    # Calculate the scale in meters per pixel.
    scale_dist_px = scale_200_px - scale_0_px
    scale_m_per_px = scale_dist_m / scale_dist_px
    print(f"Map scale: {scale_dist_m} meters / {scale_dist_px} pixels = {scale_m_per_px:.4f} m/px")

    # Calculate the horizontal distances between points in meters.
    # Side c of triangle XYZ
    dist_XY_m = math.sqrt((px_X[0] - px_Y[0])**2 + (px_X[1] - px_Y[1])**2) * scale_m_per_px
    # Side a of triangle XYZ
    dist_YZ_m = math.sqrt((px_Y[0] - px_Z[0])**2 + (px_Y[1] - px_Z[1])**2) * scale_m_per_px
    # Side b of triangle XYZ
    dist_XZ_m = math.sqrt((px_X[0] - px_Z[0])**2 + (px_X[1] - px_Z[1])**2) * scale_m_per_px
    print(f"Horizontal distance XY = {dist_XY_m:.2f} m")
    print(f"Horizontal distance YZ = {dist_YZ_m:.2f} m")
    print(f"Horizontal distance XZ = {dist_XZ_m:.2f} m\n")

    # Step 2: Determine the strike line.
    # Find point P on line YZ with the same elevation as X (120m).
    # The position of P is proportional to the elevation difference.
    print("Step 2: Finding the strike line...")
    height_ratio = (h_X - h_Y) / (h_Z - h_Y)
    dist_YP_m = dist_YZ_m * height_ratio
    print(f"A point 'P' with elevation {h_X}m exists on the line YZ.")
    print(f"The distance from Y to P is {height_ratio:.2f} of the total distance from Y to Z.")
    print(f"Calculated distance YP = {dist_YP_m:.2f} m.")
    print("The line connecting X and P is a strike line at 120m elevation.\n")

    # Step 3: Calculate the horizontal distance (run) perpendicular to the strike.
    # We calculate this by finding the altitude of triangle XYP from vertex Y to base XP.
    print("Step 3: Calculating horizontal distance (run) for dip measurement...")
    # First, find the angle at vertex Y in triangle XYZ using the Law of Cosines.
    cos_angle_Y = (dist_YZ_m**2 + dist_XY_m**2 - dist_XZ_m**2) / (2 * dist_YZ_m * dist_XY_m)
    angle_Y_rad = math.acos(cos_angle_Y)
    
    # Now, find the length of the strike line XP in triangle XYP using the Law of Cosines.
    dist_XP_m_sq = dist_XY_m**2 + dist_YP_m**2 - 2 * dist_XY_m * dist_YP_m * math.cos(angle_Y_rad)
    dist_XP_m = math.sqrt(dist_XP_m_sq)
    
    # Calculate the area of triangle XYP.
    area_XYP = 0.5 * dist_XY_m * dist_YP_m * math.sin(angle_Y_rad)
    
    # The horizontal run is the altitude of triangle XYP from vertex Y.
    # run = (2 * Area) / base
    run = (2 * area_XYP) / dist_XP_m
    print(f"The perpendicular horizontal distance from Y to the strike line XP is {run:.2f} m.\n")

    # Step 4: Calculate the dip angle.
    print("Step 4: Calculating the dip angle...")
    # The vertical rise is the elevation difference between the strike line (at h_X) and point Y.
    rise = h_X - h_Y
    print(f"Vertical rise = Elevation(X) - Elevation(Y) = {h_X} m - {h_Y} m = {rise} m")
    print(f"Horizontal run = {run:.1f} m")

    # Calculate dip angle in radians, then convert to degrees.
    dip_rad = math.atan(rise / run)
    dip_deg = math.degrees(dip_rad)
    
    print("\nThe dip angle is calculated using the formula: Dip = arctan(Rise / Run)")
    print(f"Dip = arctan({rise} m / {run:.1f} m)")
    print(f"Dip = {dip_deg:.2f} degrees")
    
    # Round to the nearest degree.
    rounded_dip = round(dip_deg)
    print(f"Rounded to the nearest degree, the dip is {rounded_dip} degrees.\n")
    
    return rounded_dip

# Execute the calculation and print the final answer.
final_answer = calculate_dip()
print(f"<<<{final_answer}>>>")

calculate_dip()