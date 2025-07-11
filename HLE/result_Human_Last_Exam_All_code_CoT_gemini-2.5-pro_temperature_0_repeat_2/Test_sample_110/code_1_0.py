import math

def solve_robot_arm_distance():
    """
    Calculates the minimum distance between the robot's finger and shoulder
    considering all geometric and collision constraints.
    """
    # Segment lengths in cm
    L1 = 40.0  # Shoulder to Elbow
    L2 = 28.0  # Elbow to Wrist
    L3 = 15.0  # Wrist to Hand
    L4 = 10.0  # Hand (finger)

    # Constraints
    clearance = 3.5
    collision_dist = 1.0

    # --- Step 1: Define minimum joint angles from clearance ---
    # Model: alpha = 2 * asin(clearance / L_shorter)
    # These are the tightest possible folds for the wrist and finger.
    alpha_w_min = 2 * math.asin(clearance / L3)  # L3 is shorter than L2
    alpha_f_min = 2 * math.asin(clearance / L4)  # L4 is shorter than L3

    # --- Step 2: Find the elbow angle that avoids collision ---
    # The problem is to find alpha_e that satisfies the collision constraint.
    # We need to solve for alpha_e in the equation where the lowest point of
    # segment L3 is exactly `collision_dist` above the axis of segment L1.
    # The equation for the y-coordinate of the hand joint (P3) is:
    # y_P3 = L2*sin(a_e) - L3*sin(a_e + a_w_min)
    # We set y_P3 = collision_dist and solve for a_e.
    # This is a transcendental equation of the form A*sin(x) + B*cos(x) = C.
    # Let's solve it numerically using a simple search.
    
    alpha_e = 0.0  # Initial guess
    # Search for the correct alpha_e
    for i in range(5000):
        ae = i * 0.0001  # Increment alpha_e
        # y_P3 = L2*sin(ae) - L3*sin(ae + alpha_w_min)
        # We need to expand sin(ae + aw)
        y_p3 = L2 * math.sin(ae) - L3 * (math.sin(ae) * math.cos(alpha_w_min) + math.cos(ae) * math.sin(alpha_w_min))
        if y_p3 >= collision_dist:
            alpha_e = ae
            break

    # Use the found alpha_e, and the minimums for the other joints
    alpha_w = alpha_w_min
    alpha_f = alpha_f_min

    # --- Step 3: Calculate final tip coordinates with the valid angles ---
    # We use forward kinematics with the corrected elbow angle.
    # The arm folds back, so angles are subtracted from pi.
    # Absolute angle of L2
    angle2 = math.pi - alpha_e
    # Absolute angle of L3
    angle3 = angle2 + (math.pi - alpha_w)
    # Absolute angle of L4
    angle4 = angle3 + (math.pi - alpha_f)

    # P1 (Elbow)
    p1_x = L1
    p1_y = 0

    # P2 (Wrist)
    p2_x = p1_x + L2 * math.cos(angle2)
    p2_y = p1_y + L2 * math.sin(angle2)

    # P3 (Hand joint)
    p3_x = p2_x + L3 * math.cos(angle3)
    p3_y = p2_y + L3 * math.sin(angle3)

    # P_tip (Finger tip)
    p_tip_x = p3_x + L4 * math.cos(angle4)
    p_tip_y = p3_y + L4 * math.sin(angle4)

    # --- Step 4: Calculate and print the final distance ---
    distance = math.sqrt(p_tip_x**2 + p_tip_y**2)
    
    print("The final position of the finger tip is ({:.2f} cm, {:.2f} cm).".format(p_tip_x, p_tip_y))
    print("The distance from the shoulder is calculated by the equation:")
    print("distance = sqrt({:.2f}^2 + {:.2f}^2) = {:.2f} cm".format(p_tip_x, p_tip_y, distance))

solve_robot_arm_distance()