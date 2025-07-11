import math

def calculate_min_distance():
    """
    Calculates the minimum distance from the robot's finger to its shoulder
    based on geometric constraints.
    """

    # Segment lengths in cm
    l_se = 40  # Shoulder to Elbow
    l_ew = 28  # Elbow to Wrist
    l_wh = 15  # Wrist to Hand
    l_hf = 10  # Hand to Finger

    # Step 1: Determine the minimum internal joint angles from constraints.
    # The minimum distance between non-adjacent segments is 1 cm.
    # We model this as: separation_distance = separator_length * sin(internal_angle)
    
    # Constraint: dist(SE, WH) >= 1. Separator is EW.
    # 28 * sin(alpha_E) = 1
    # We check if the argument to asin is valid.
    if 1 / l_ew > 1:
        print("Error: Impossible geometry for elbow constraint.")
        return

    alpha_e = math.asin(1 / l_ew)  # Elbow internal angle in radians

    # Constraint: dist(EW, HF) >= 1. Separator is WH.
    # 15 * sin(alpha_W) = 1
    if 1 / l_wh > 1:
        print("Error: Impossible geometry for wrist constraint.")
        return
    alpha_w = math.asin(1 / l_wh)  # Wrist internal angle in radians
    
    # For the hand joint, there's no non-adjacent segment constraint.
    # To get the finger closest to the shoulder, we assume it folds completely.
    alpha_h = 0.0  # Hand internal angle in radians

    # Step 2: Calculate the final coordinates of the fingertip.
    # Start with the shoulder at (0,0) and the elbow at (L_SE, 0).
    x, y = l_se, 0.0
    
    # Start with the first segment's angle (global frame)
    current_angle = 0.0
    
    # Add Elbow-Wrist segment
    # The joint folds, so the segment turns by (pi - alpha) relative to the previous one.
    current_angle += math.pi - alpha_e
    x += l_ew * math.cos(current_angle)
    y += l_ew * math.sin(current_angle)
    
    # Add Wrist-Hand segment
    current_angle += math.pi - alpha_w
    x += l_wh * math.cos(current_angle)
    y += l_wh * math.sin(current_angle)

    # Add Hand-Finger segment
    current_angle += math.pi - alpha_h
    x += l_hf * math.cos(current_angle)
    y += l_hf * math.sin(current_angle)

    # Step 3: Calculate the distance from the origin (shoulder)
    distance = math.sqrt(x**2 + y**2)

    # The equation for the final distance can be written out for clarity
    print("The final position of the fingertip is approximately ({:.2f}, {:.2f}) cm.".format(x, y))
    print("The minimum distance is the hypotenuse of this position:")
    print("sqrt({:.2f}^2 + {:.2f}^2) = {:.2f} cm".format(x, y, distance))


calculate_min_distance()