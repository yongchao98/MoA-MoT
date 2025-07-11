import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron
    from an isosceles right triangle based on Synge's formula.
    """

    # Step 1: Define the angles of the isosceles right triangle in radians.
    # The angle opposite the hypotenuse is the right angle.
    angle_opposite_hypotenuse = math.pi / 2
    # The two angles opposite the legs are equal.
    angle_opposite_leg = math.pi / 4

    # Step 2: Use Synge's formula P_i = (pi - A_i) / (2*pi) to find the probabilities.
    # A_i is the angle opposite side i.

    # Probability of escaping through the hypotenuse (H).
    prob_hypotenuse = (math.pi - angle_opposite_hypotenuse) / (2 * math.pi)

    # Probability of escaping through one leg (L1 or L2).
    prob_one_leg = (math.pi - angle_opposite_leg) / (2 * math.pi)

    # The total probability of escaping through either of the two legs.
    prob_legs_total = 2 * prob_one_leg

    # Step 3: Calculate the difference.
    difference = prob_hypotenuse - prob_legs_total

    # Step 4: Output the results and the final equation.
    print("This script calculates the required probability difference using an exact analytical formula.")
    print("-" * 70)
    print(f"Angle opposite the hypotenuse: {math.degrees(angle_opposite_hypotenuse):.0f} degrees (pi/2 radians)")
    print(f"Angle opposite each leg: {math.degrees(angle_opposite_leg):.0f} degrees (pi/4 radians)")
    print("-" * 70)

    print("Calculating the probability of escaping through the hypotenuse, P(H):")
    print(f"P(H) = (pi - pi/2) / (2*pi) = {prob_hypotenuse:.4f}")

    print("\nCalculating the total probability of escaping through the two legs, P(L_total):")
    print(f"P(L_total) = 2 * (pi - pi/4) / (2*pi) = {prob_legs_total:.4f}")

    print("\nCalculating the final difference P(H) - P(L_total):")
    print(f"{prob_hypotenuse:.4f} - {prob_legs_total:.4f} = {difference:.4f}")


solve_escape_probability()
<<< -0.5 >>>