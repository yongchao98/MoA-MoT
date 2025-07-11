import math

def calculate_final_distance():
    """
    Calculates the minimum distance from a robot arm's finger to its shoulder.
    """
    # Segment lengths in cm
    l1 = 40.0  # Shoulder to Elbow
    l2 = 28.0  # Elbow to Wrist
    l3 = 15.0  # Wrist to Hand
    l4 = 10.0  # Hand to Finger
    
    # Clearance constraint in cm
    c = 3.5

    # --- Step 1: Calculate the minimum folding angle for each joint ---
    # Model: l_moving * sin(alpha) = c
    # We need to handle potential math domain errors if c > l
    if c > l2 or c > l3 or c > l4:
        raise ValueError("Clearance is too large to be physically possible for a segment.")

    # Elbow joint (angle between l1 and l2)
    sin_alpha_elbow = c / l2
    cos_alpha_elbow = math.sqrt(1 - sin_alpha_elbow**2)
    alpha_elbow = math.asin(sin_alpha_elbow)

    # Wrist joint (angle between l2 and l3)
    sin_alpha_wrist = c / l3
    cos_alpha_wrist = math.sqrt(1 - sin_alpha_wrist**2)
    alpha_wrist = math.asin(sin_alpha_wrist)

    # Hand joint (angle between l3 and l4)
    sin_alpha_hand = c / l4
    cos_alpha_hand = math.sqrt(1 - sin_alpha_hand**2)
    alpha_hand = math.asin(sin_alpha_hand)

    # --- Step 2: Calculate coordinates using a 2D kinematic chain (zig-zag fold) ---
    # Shoulder is at the origin (0, 0)
    # P_s = (0, 0)
    
    # Elbow position
    p_elbow_x = l1
    p_elbow_y = 0.0

    # Wrist position (folds back from elbow)
    # Global angle of segment 2 is pi - alpha_elbow
    p_wrist_x = p_elbow_x - l2 * cos_alpha_elbow
    p_wrist_y = p_elbow_y + l2 * sin_alpha_elbow
    
    # Hand position (folds forward from wrist)
    # Global angle of segment 3 is (pi - alpha_elbow) - (pi - alpha_wrist) = alpha_wrist - alpha_elbow
    # We use trigonometric identities for cos(A-B) and sin(A-B)
    cos_aw_minus_ae = cos_alpha_wrist * cos_alpha_elbow + sin_alpha_wrist * sin_alpha_elbow
    sin_aw_minus_ae = sin_alpha_wrist * cos_alpha_elbow - cos_alpha_wrist * sin_alpha_elbow
    p_hand_x = p_wrist_x + l3 * cos_aw_minus_ae
    p_hand_y = p_wrist_y + l3 * sin_aw_minus_ae
    
    # Finger position (folds back from hand)
    # Global angle is (alpha_wrist - alpha_elbow) + (pi - alpha_hand)
    # We use identities for cos(pi+X) = -cos(X) and sin(pi+X) = -sin(X) where X = aw-ae-ah
    # First, find cos(aw-ae-ah) and sin(aw-ae-ah)
    cos_aw_ae = cos_aw_minus_ae
    sin_aw_ae = sin_aw_minus_ae
    cos_X = cos_aw_ae * cos_alpha_hand + sin_aw_ae * sin_alpha_hand
    sin_X = sin_aw_ae * cos_alpha_hand - cos_aw_ae * sin_alpha_hand
    
    # Add the vector for the final segment
    p_finger_x = p_hand_x - l4 * cos_X
    p_finger_y = p_hand_y - l4 * sin_X

    # --- Step 3: Calculate the final distance from the shoulder (origin) ---
    distance = math.sqrt(p_finger_x**2 + p_finger_y**2)
    
    print("This calculation determines the final position of the robot's finger after folding as much as possible based on a clearance constraint.")
    print("The arm segments have lengths L1={}, L2={}, L3={}, L4={} cm.".format(l1, l2, l3, l4))
    print("The folding at each joint is limited by a clearance value of c={} cm.".format(c))
    print("The final position of the finger is ({:.2f}, {:.2f}) cm relative to the shoulder.".format(p_finger_x, p_finger_y))
    print("The distance is the hypotenuse of these coordinates:")
    print("Distance = sqrt({:.2f}^2 + {:.2f}^2) = {:.2f} cm".format(p_finger_x, p_finger_y, distance))

calculate_final_distance()