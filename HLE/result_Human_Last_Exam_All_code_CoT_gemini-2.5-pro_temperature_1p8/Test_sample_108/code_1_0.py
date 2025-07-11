import math

def solve_robot_arm_distance():
    """
    Calculates the minimum distance from the robot arm's finger to its shoulder.
    """
    # Segment lengths in cm
    L1 = 40  # Shoulder to Elbow
    L2 = 28  # Elbow to Wrist
    L3 = 15  # Wrist to Hand
    L4 = 10  # Hand (Finger)

    # Collision constraint: minimum distance between non-adjacent segments
    min_dist = 1.0

    # Based on geometric analysis, the minimum required vertical separation 'D'
    # is dictated by the collision constraint between segments L2 and L4.
    # We found that to keep the distance between L2 and L4 >= 1 cm,
    # the separation D must be at least 28/5 cm.
    # d(L2,L4) = 5*D/28 >= 1  => D >= 28/5
    D = 28.0 / 5.0

    # Let's verify this satisfies the other main constraint, d(L2, a point on L3), which was 15*D/28 >= 1
    # 15 * (28/5) / 28 = 15/5 = 3. This is >= 1, so the D=28/5 is the correct minimal separation.

    print(f"Segment lengths (cm): L1={L1}, L2={L2}, L3={L3}, L4={L4}")
    print(f"Minimum required separation (D) to avoid collision: {D:.4f} cm")

    # Now, calculate the coordinates of the fingertip (xF, yF)
    # The configuration is: L1(+x), L2(-x), L3(-x), L4(+x), with segments on y=0 or y=D.

    # y-coordinate of the fingertip is the separation D
    yF = D
    print(f"Final finger y-coordinate: {yF:.4f} cm")


    # x-coordinate calculation:
    # x-coordinate of the wrist (W)
    # xW = L1 - sqrt(L2^2 - D^2)
    xW = L1 - math.sqrt(L2**2 - D**2)
    
    # x-coordinate of the hand (H)
    # xH = xW - L3
    xH = xW - L3

    # x-coordinate of the finger (F)
    # xF = xH + L4
    xF = xH + L4
    
    print(f"Final finger x-coordinate: {xF:.4f} cm")

    # Final distance is the hypotenuse of a right triangle with sides xF and yF
    final_distance = math.sqrt(xF**2 + yF**2)

    # Output the final equation with numbers
    print(f"\nFinal Distance = sqrt( ( {L1} - sqrt({L2}^2 - {D:.4f}^2) - {L3} + {L4} )^2 + {D:.4f}^2 )")
    print(f"Final Distance = sqrt( {xF:.4f}^2 + {yF:.4f}^2 )")
    print(f"Final Distance = sqrt( {xF**2:.4f} + {yF**2:.4f} )")
    print(f"Final Distance = sqrt( {xF**2 + yF**2:.4f} )")
    print(f"\nCalculated minimum distance from finger to shoulder: {final_distance:.4f} cm")

solve_robot_arm_distance()
<<<A>>>