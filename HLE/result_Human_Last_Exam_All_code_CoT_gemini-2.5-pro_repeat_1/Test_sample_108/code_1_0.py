import math

def calculate_min_distance():
    """
    Calculates the minimum distance of a robot arm's finger from its shoulder.
    """
    # --- Step 1: Define constants ---
    L1 = 40.0  # Shoulder to Elbow length
    L2 = 28.0  # Elbow to Wrist length
    L3 = 15.0  # Wrist to Hand length
    L4 = 10.0  # Hand to Finger length
    clearance = 3.5
    
    print("Arm Segment Lengths:")
    print(f"L1 (Shoulder-Elbow): {L1} cm")
    print(f"L2 (Elbow-Wrist): {L2} cm")
    print(f"L3 (Wrist-Hand): {L3} cm")
    print(f"L4 (Hand-Finger): {L4} cm")
    print("-" * 20)

    # --- Step 2: Calculate joint angle deviations ---
    # The minimum deviation angle 'delta' from a 180-degree fold is calculated
    # based on the clearance constraint: sin(delta) = clearance / next_segment_length
    delta_elbow = math.asin(clearance / L2)
    delta_wrist = math.asin(clearance / L3)
    delta_hand = math.asin(clearance / L4)

    # The joint angles are pi - delta, representing the maximum fold.
    angle_elbow_joint = math.pi - delta_elbow
    angle_wrist_joint = math.pi - delta_wrist
    angle_hand_joint = math.pi - delta_hand

    # --- Step 3: Calculate vector for each segment in a spiral configuration ---
    # We start with L1 along the positive x-axis.
    # We then sequentially add the other vectors, each rotated relative to the previous one.
    
    # Position starts at the origin (shoulder)
    final_x = 0.0
    final_y = 0.0
    current_angle = 0.0

    # Vector for L1
    v1_x = L1 * math.cos(current_angle)
    v1_y = L1 * math.sin(current_angle)
    final_x += v1_x
    final_y += v1_y

    # Vector for L2
    current_angle += angle_elbow_joint
    v2_x = L2 * math.cos(current_angle)
    v2_y = L2 * math.sin(current_angle)
    final_x += v2_x
    final_y += v2_y

    # Vector for L3
    current_angle += angle_wrist_joint
    v3_x = L3 * math.cos(current_angle)
    v3_y = L3 * math.sin(current_angle)
    final_x += v3_x
    final_y += v3_y

    # Vector for L4
    current_angle += angle_hand_joint
    v4_x = L4 * math.cos(current_angle)
    v4_y = L4 * math.sin(current_angle)
    final_x += v4_x
    final_y += v4_y

    # --- Step 4: Calculate the final distance ---
    distance = math.sqrt(final_x**2 + final_y**2)

    # --- Step 5: Print the final equation and result ---
    print("Final Position Calculation (Sum of Vectors):")
    print(f"X = {v1_x:.2f} + ({v2_x:.2f}) + ({v3_x:.2f}) + ({v4_x:.2f}) = {final_x:.2f} cm")
    print(f"Y = {v1_y:.2f} + {v2_y:.2f} + ({v3_y:.2f}) + {v4_y:.2f} = {final_y:.2f} cm")
    print("-" * 20)
    print("Final Distance Calculation (Magnitude of final vector):")
    print(f"Distance = sqrt( ({final_x:.2f})^2 + ({final_y:.2f})^2 )")
    print(f"Distance = {distance:.2f} cm")

# Run the calculation
calculate_min_distance()