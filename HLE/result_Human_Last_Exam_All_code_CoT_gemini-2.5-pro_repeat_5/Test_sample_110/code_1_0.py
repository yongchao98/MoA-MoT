import math

def calculate_robot_arm_distance():
    """
    Calculates the minimum distance between the robot's finger and shoulder
    based on a simplified perpendicular folding model.
    """
    # Segment lengths
    L1 = 40  # Shoulder to Elbow
    L2 = 28  # Elbow to Wrist
    L3 = 15  # Wrist to Hand
    L4 = 10  # Hand to Finger

    # Step 1: Calculate the effective length from shoulder to wrist (d_sw)
    # This assumes the first two segments fold back on each other.
    d_sw = abs(L1 - L2)

    # Step 2: Calculate the effective length from wrist to finger (d_wf)
    # This assumes the last two segments fold back on each other.
    d_wf = abs(L3 - L4)

    # Step 3: Combine the effective lengths.
    # We assume the two folded sections (S-W and W-F) arrange themselves
    # perpendicularly to avoid self-collision, which is a complex constraint.
    # The final distance is the hypotenuse of a right triangle with sides d_sw and d_wf.
    final_distance = math.sqrt(d_sw**2 + d_wf**2)
    
    print("This calculation models the arm folding into two main sections.")
    print("The first section is from the shoulder to the wrist.")
    print(f"Its ideal folded length is |L1 - L2| = |{L1} - {L2}| = {d_sw} cm.")
    print("\nThe second section is from the wrist to the finger.")
    print(f"Its ideal folded length is |L3 - L4| = |{L3} - {10}| = {d_wf} cm.")
    print("\nAssuming these two sections are perpendicular to avoid collision,")
    print("the final distance is calculated using the Pythagorean theorem:")
    print(f"Distance = sqrt({d_sw}^2 + {d_wf}^2)")
    print(f"Distance = sqrt({d_sw**2} + {d_wf**2})")
    print(f"Distance = sqrt({d_sw**2 + d_wf**2})")
    print(f"Final calculated distance: ~{final_distance:.2f} cm")

calculate_robot_arm_distance()