import math

def solve_conic_asymptote_angles(alpha_deg, beta_deg, gamma_deg, delta_deg):
    """
    Calculates the angles between the asymptotes of a specific conic and the line BC.

    Args:
        alpha_deg: Angle alpha of triangle ABC in degrees.
        beta_deg: Angle beta of triangle ABC in degrees.
        gamma_deg: Angle gamma of triangle ABC in degrees.
        delta_deg: Angle between line BC and line l in degrees.

    Returns:
        A tuple of two angles in radians.
    """
    # Convert delta from degrees to radians for calculation
    delta = math.radians(delta_deg)

    # The angles of the asymptotes with the line BC depend only on delta.
    # The first angle is -delta.
    angle1 = -delta
    
    # The conic is a rectangular hyperbola, so its asymptotes are perpendicular.
    # The second angle is -delta + pi/2.
    angle2 = -delta + math.pi / 2

    # We can normalize the angles to be in the range [0, pi) as they represent lines
    # The set of angles is {(-delta) mod pi, (-delta + pi/2) mod pi}
    
    print("This problem is a highly complex geometry question whose solution is surprisingly simple.")
    print("The angles of the asymptotes with the line BC are independent of the triangle's angles (alpha, beta, gamma).")
    print("They only depend on the angle delta between line BC and line l.")
    print("\nThe formulas for the angles (in radians) are:")
    print("angle1 = -delta")
    print("angle2 = pi/2 - delta")
    
    print(f"\nGiven delta = {delta_deg} degrees ({delta:.4f} radians):")
    print(f"The first angle is -delta.")
    print(f"Equation: angle1 = -{delta:.4f}")
    print(f"Result: {angle1:.4f} radians, which is {math.degrees(angle1):.2f} degrees.")
    
    print(f"\nThe second angle is pi/2 - delta.")
    print(f"Equation: angle2 = {math.pi:.4f} / 2 - {delta:.4f}")
    print(f"Result: {angle2:.4f} radians, which is {math.degrees(angle2):.2f} degrees.")

# Example usage with sample values for the angles
# Let's use a standard triangle and a non-trivial angle for delta.
alpha_sample = 60  # degrees
beta_sample = 70   # degrees
gamma_sample = 50   # degrees
delta_sample = 30  # degrees

solve_conic_asymptote_angles(alpha_sample, beta_sample, gamma_sample, delta_sample)
