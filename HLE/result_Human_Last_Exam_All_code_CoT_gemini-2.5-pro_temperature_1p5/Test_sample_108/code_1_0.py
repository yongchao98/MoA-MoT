import math

def solve_robot_arm_distance():
    """
    Calculates the minimum distance from the robot's shoulder to its finger.
    """
    l1 = 40.0  # Shoulder to elbow
    l2 = 28.0  # Elbow to wrist
    l3 = 15.0  # Wrist to hand
    l4 = 10.0  # Hand to finger
    clearance = 1.0 # Minimum distance between non-adjacent segments in cm

    # We model a 2D folding configuration where L1 is on the x-axis, and the
    # other segments fold back on a parallel line at a distance 'd'.
    # This configuration is L1 ->, L2 <-, L3 <-, L4 ->.
    # The ideal, co-linear distance would be 40 - 28 - 15 + 10 = 7 cm.
    # Now, we account for the 1cm clearance for non-adjacent segments.

    # The most restrictive constraint is the distance between L2 and L4.
    # Let's find the minimum separation 'd' required.
    # The minimum distance between segment L2 and segment L4 in this parallel
    # model can be shown to be d * |L3 - L4| / L2.
    # So, we need: d * |L3 - L4| / L2 >= clearance
    # d >= clearance * L2 / |L3 - L4|

    min_d = clearance * l2 / abs(l3 - l4)

    # We also need d >= clearance for the L1 vs L3 and L1 vs L4 constraints.
    # Since L2 / |L3 - L4| = 28 / 5 = 5.6, this is the dominant constraint.
    required_d = max(clearance, min_d)

    # Now calculate the final fingertip position (P4) with this separation d.
    # P4.x = L1_x - L2_projection_x - L3_x + L4_x
    # The projection of L2 onto the x-axis is sqrt(L2^2 - d^2)
    p4_x = l1 - math.sqrt(l2**2 - required_d**2) - l3 + l4
    p4_y = required_d

    # The final distance is the magnitude of the vector to P4.
    final_distance = math.sqrt(p4_x**2 + p4_y**2)

    print(f"Segment lengths (cm): L1={l1}, L2={l2}, L3={l3}, L4={l4}")
    print(f"Required clearance between non-adjacent segments: {clearance} cm")
    print("\nAnalysis of the 'L1-L2-L3+L4' folding configuration:")
    print(f"To prevent collision between segment L2 and L4, the parallel separation 'd' must be at least: {required_d:.2f} cm")
    print(f"With this separation, the fingertip is at coordinate (x,y) relative to the shoulder.")
    print(f"P4_x = L1 - sqrt(L2^2 - d^2) - L3 + L4 = {l1} - sqrt({l2}^2 - {required_d:.2f}^2) - {l3} + {l4} = {p4_x:.2f} cm")
    print(f"P4_y = d = {required_d:.2f} cm")
    print("\nThe final distance from shoulder to finger is the length of the vector to P4 (x,y):")
    print(f"Distance = sqrt(P4_x^2 + P4_y^2) = sqrt({p4_x:.2f}^2 + {required_d:.2f}^2) = {final_distance:.2f} cm")

solve_robot_arm_distance()