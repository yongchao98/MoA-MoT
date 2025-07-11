import numpy as np
from scipy.optimize import brentq

# Define the lengths of the robot arm segments
L1 = 40.0  # Shoulder to Elbow
L2 = 28.0  # Elbow to Wrist
L3 = 15.0  # Wrist to Hand
L4 = 10.0  # Hand (finger)

# Define the collision constraint
MIN_CLEARANCE = 1.0

# This function represents the core geometric constraint of the problem.
# It relates the distance from shoulder to wrist (d) to the geometry of the
# L1-L2 triangle and the collision clearance constraint.
# We want to find the root of this function.
def constraint_equation(d):
    """
    Calculates the value of the constraint function for a given distance 'd'.
    The function is f(d) = x^2 + y^2 - d^2, where x and y are derived from
    the arm's geometry and constraints. The root of this function is the
    distance 'd' that satisfies all conditions simultaneously.

    Args:
        d (float): The distance from the shoulder (P0) to the wrist (P2).

    Returns:
        float: The value of the constraint function. Should be 0 for a valid configuration.
    """
    # d must be greater than L3 for this model to be valid
    if d <= L3:
        # Return a large number to guide the solver away from this region.
        # This occurs if the hand P3 would be at or beyond the shoulder P0.
        return 1e9

    # From the law of cosines in triangle P0-P1-P2, or by solving the circle equation,
    # we can find the x-coordinate of the wrist P2.
    # (x-L1)^2 + y^2 = L2^2  and x^2 + y^2 = d^2
    # d^2 - 2*L1*x + L1^2 = L2^2 => x = (d^2 + L1^2 - L2^2) / (2*L1)
    # Using the provided numbers: x = (d^2 + 40^2 - 28^2) / (2*40) = (d^2 + 816) / 80
    x = (d**2 + L1**2 - L2**2) / (2 * L1)

    # The y-coordinate of the hand joint P3 must be >= MIN_CLEARANCE.
    # The most restrictive case is when this constraint is active (equal to clearance).
    # The y-coord of P3 is y_p3 = y_p2 * (1 - L3/d).
    # So, MIN_CLEARANCE = y_p2 * (1 - L3/d).
    # This gives y_p2 = MIN_CLEARANCE / (1 - L3/d) = MIN_CLEARANCE * d / (d - L3)
    y = MIN_CLEARANCE * d / (d - L3)

    # The final constraint equation is from the definition of d: d^2 = x^2 + y^2
    # So, f(d) = x^2 + y^2 - d^2 = 0
    return x**2 + y**2 - d**2

# We need to find the root of the constraint_equation.
# From manual inspection, we know the solution for d is between 15 and 20.
# Let's use a root-finding algorithm to find the precise value.
try:
    # Find the minimum valid distance 'd' from shoulder to wrist
    d_02_min = brentq(constraint_equation, L3 + 1e-6, L1 - L2 + L3)
    
    # The final distance from shoulder to finger is d_02_min - |L3 - L4|
    dist_finger_to_shoulder = d_02_min - abs(L3 - L4)
    
    print("This solution is derived by modeling the robot arm's geometry and applying the non-adjacent segment collision constraint.")
    print("The key constraints are:")
    print(f"1. The length of the shoulder-to-elbow segment (L1): {L1} cm")
    print(f"2. The length of the elbow-to-wrist segment (L2): {L2} cm")
    print(f"3. The length of the wrist-to-hand segment (L3): {L3} cm")
    print(f"4. The length of the hand-to-finger segment (L4): {L4} cm")
    print(f"5. The minimum clearance between non-adjacent segments: {MIN_CLEARANCE} cm")
    print("\nTo minimize the final distance, we fold the arm back on itself as much as the constraints allow.")
    print("The most restrictive constraint is the 1 cm clearance between the first segment (L1) and the third segment (L3).")
    print("\nLet d be the distance from the shoulder to the wrist.")
    print("Solving the geometric constraint equations numerically gives the minimum required value for d.")
    print(f"\nMinimum required distance from shoulder to wrist (d): {d_02_min:.4f} cm")
    print("The final distance to the finger is calculated by folding the last two segments back.")
    print(f"Final Distance = d - |L3 - L4|")
    print(f"Final Distance = {d_02_min:.4f} - |{L3} - {L4}|")
    print(f"Final Distance = {d_02_min:.4f} - {abs(L3 - L4)}")
    print(f"Final Distance = {dist_finger_to_shoulder:.4f} cm")
    print(f"\nThe closest answer is {dist_finger_to_shoulder:.2f} cm.")

except ValueError:
    print("Could not find a solution in the given interval.")
