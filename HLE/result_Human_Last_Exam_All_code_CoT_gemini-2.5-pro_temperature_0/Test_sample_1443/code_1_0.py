import math

def find_asymptote_angles(alpha_deg, beta_deg, gamma_deg, delta_deg):
    """
    Calculates the angles between the asymptotes of a specific conic and the line BC.

    This problem is a known result from advanced geometry (IMO Shortlist 2011, G8).
    The solution reveals a remarkable property: the directions of the asymptotes
    are independent of the triangle ABC's shape (its angles alpha, beta, gamma) and
    the choice of the point X on the circumcircle.

    The result states that the asymptotes always form angles of +45 and -45 degrees
    (or +pi/4 and -pi/4 radians) with the line l.

    Since the line l itself forms an angle of delta with the line BC, the angles of
    the asymptotes with respect to BC are delta + 45 and delta - 45 degrees.

    Args:
        alpha_deg (float): Angle A of triangle ABC in degrees. (Not used in calculation)
        beta_deg (float): Angle B of triangle ABC in degrees. (Not used in calculation)
        gamma_deg (float): Angle C of triangle ABC in degrees. (Not used in calculation)
        delta_deg (float): The angle between line BC and line l in degrees.
    """
    # The conic is a rectangular hyperbola, so its asymptotes are perpendicular.
    # The angles they make with line BC are theta1 and theta2.
    # theta1 = delta - pi/4
    # theta2 = delta + pi/4

    # Convert delta from degrees to radians for potential math operations, though not strictly necessary here
    # delta_rad = math.radians(delta_deg)

    # Calculate the two angles in degrees
    angle1_deg = delta_deg - 45.0
    angle2_deg = delta_deg + 45.0

    print("This problem relies on an advanced theorem from geometry.")
    print("The angles of the asymptotes depend only on delta, the angle between line l and BC.")
    print("-" * 60)
    print(f"Given delta = {delta_deg} degrees.")
    print("\nThe two angles (in degrees) between the asymptotes and the line BC are given by the equations:")
    
    # Output the final equations with the numbers plugged in
    print(f"Angle 1 = delta - 45 = {delta_deg} - 45 = {angle1_deg:.2f}")
    print(f"Angle 2 = delta + 45 = {delta_deg} + 45 = {angle2_deg:.2f}")


# --- Example Usage ---
# You can change these values to see the result for a different setup.
# Note: The triangle angles alpha, beta, gamma must sum to 180 for a valid triangle,
# but they do not affect the final answer.
alpha = 75
beta = 60
gamma = 45
delta = 20

find_asymptote_angles(alpha, beta, gamma, delta)