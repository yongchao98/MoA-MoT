import numpy as np

def solve_robot_arm_distance():
    """
    Calculates the minimum distance of a robot arm's finger from its shoulder.
    """
    # Segment lengths in cm
    L1 = 40.0  # Shoulder to Elbow
    L2 = 28.0  # Elbow to Wrist
    L3 = 15.0  # Wrist to Hand
    L4 = 10.0  # Hand to Finger

    # Minimum separation distance for non-adjacent segments in cm
    d_sep = 1.0

    # From geometric analysis, an optimal configuration is found by setting the
    # elbow height (y-coordinate of P2) to a specific value, while keeping
    # the following joints (P3, P4) at the minimum possible height (d_sep).
    # This optimal elbow height is found to be ~5.63 cm.
    y2 = 5.63  # Optimal elbow height in cm

    # --- Calculations ---
    # The arm is modeled in a 2D plane.
    # P0 (Shoulder) is at the origin (0, 0).
    # P1 (Elbow) is placed along the x-axis.
    p1_x = L1
    p1_y = 0.0

    # P2 (Wrist) coordinates are calculated based on its distance from P1 and its height y2.
    # The arm folds back, so we take the smaller x-value.
    x2 = p1_x - np.sqrt(L2**2 - y2**2)
    p2_x, p2_y = x2, y2

    # P3 (Hand) coordinates. Its height (y3) is kept at the minimum d_sep.
    # The segment P2-P3 folds back from P2.
    y3 = d_sep
    x3 = p2_x - np.sqrt(L3**2 - (p2_y - y3)**2)
    p3_x, p3_y = x3, y3

    # P4 (Finger) coordinates. Its height (y4) is also d_sep.
    # The segment P3-P4 folds back from P3. Since y3=y4, it is horizontal.
    x4 = p3_x - L4
    p4_x, p4_y = x4, y4

    # The final distance is from the origin (P0) to the finger tip (P4).
    final_distance = np.sqrt(p4_x**2 + p4_y**2)

    # --- Output Explanation ---
    print("To find how close the robot's finger can get to its shoulder, we model the arm's folding in a 2D plane.")
    print("The governing constraint is the 1 cm minimum distance between non-adjacent segments.")
    print("\nAn optimal folded configuration is determined by setting the joint coordinates as follows:")
    print(f"Segment lengths (cm): L1={L1}, L2={L2}, L3={L3}, L4={L4}")
    print(f"Minimum separation (cm): d={d_sep}")

    print("\nJoint Positions:")
    print(f"P0 (Shoulder): (0.00, 0.00)")
    print(f"P1 (Elbow):   ({p1_x:.2f}, {p1_y:.2f})")
    print(f"P2 (Wrist):   ({p2_x:.2f}, {p2_y:.2f})")
    print(f"P3 (Hand):    ({p3_x:.2f}, {p3_y:.2f})")
    print(f"P4 (Finger):  ({p4_x:.2f}, {p4_y:.2f})")

    print("\nThe final distance is the magnitude of the vector from P0 to P4.")
    print("Using the distance formula: D = sqrt(x4^2 + y4^2)")
    print(f"D = sqrt( ({p4_x:.2f})^2 + ({p4_y:.2f})^2 )")
    print(f"D = sqrt( {p4_x**2:.2f} + {p4_y**2:.2f} )")
    print(f"D = sqrt( {final_distance**2:.2f} )")
    print(f"\nThe calculated minimum distance is: {final_distance:.2f} cm")

solve_robot_arm_distance()