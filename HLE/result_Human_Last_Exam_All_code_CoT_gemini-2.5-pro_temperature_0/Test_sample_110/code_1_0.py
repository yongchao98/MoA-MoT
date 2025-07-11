import math

def solve_robot_arm_distance():
    """
    Calculates the distance from the robot's finger to its shoulder
    based on a simplified planar folding model.
    """
    # Segment lengths
    l1 = 40  # Shoulder to Elbow
    l2 = 28  # Elbow to Wrist
    l3 = 15  # Wrist to Hand
    l4 = 10  # Hand to Finger

    # Clearance constraint
    clearance = 1.0

    # --- Step-by-step construction of the folded arm's joint coordinates ---

    # P0: Shoulder at the origin
    p0_x, p0_y = 0.0, 0.0
    print(f"Shoulder (P0) is at: ({p0_x:.2f}, {p0_y:.2f}) cm")

    # P1: Elbow position. L1 is along the x-axis.
    p1_x, p1_y = l1, 0.0
    print(f"Elbow (P1) is at: ({p1_x:.2f}, {p1_y:.2f}) cm")

    # P2: Wrist position.
    # L2 folds back from L1. To keep L3 (starting at P2) 1cm away from L1 (at y=0),
    # P2 must have a y-coordinate of at least 1cm. We assume it's exactly 1cm.
    # This creates a right triangle with hypotenuse L2 and vertical side 'clearance'.
    dy_l2 = clearance
    dx_l2 = math.sqrt(l2**2 - dy_l2**2)
    p2_x = p1_x - dx_l2
    p2_y = p1_y + dy_l2
    print(f"Wrist (P2) is at: ({p2_x:.2f}, {p2_y:.2f}) cm")

    # P3: Hand joint position.
    # L3 folds forward from L2. To ensure clearance between L2 and L4,
    # we model L3 as also providing a vertical lift of 1cm.
    dy_l3 = clearance
    dx_l3 = math.sqrt(l3**2 - dy_l3**2)
    p3_x = p2_x + dx_l3
    p3_y = p2_y + dy_l3
    print(f"Hand (P3) is at: ({p3_x:.2f}, {p3_y:.2f}) cm")

    # P4: Fingertip position.
    # L4 folds back from P3. For this simplified model, we assume it folds
    # horizontally to maximize its backward progress.
    p4_x = p3_x - l4
    p4_y = p3_y
    print(f"Finger (P4) is at: ({p4_x:.2f}, {p4_y:.2f}) cm")

    # Calculate the final distance from the origin (shoulder) to the fingertip
    final_distance = math.sqrt(p4_x**2 + p4_y**2)

    # Print the final equation
    print("\nFinal distance calculation:")
    print(f"Distance = sqrt( (P4_x - P0_x)^2 + (P4_y - P0_y)^2 )")
    print(f"Distance = sqrt( ({p4_x:.2f} - {p0_x:.2f})^2 + ({p4_y:.2f} - {p0_y:.2f})^2 )")
    print(f"Distance = sqrt( {p4_x**2:.2f} + {p4_y**2:.2f} )")
    print(f"Distance = sqrt( {p4_x**2 + p4_y**2:.2f} )")
    print(f"\nThe final calculated distance is: {final_distance:.2f} cm")

solve_robot_arm_distance()
<<<F>>>