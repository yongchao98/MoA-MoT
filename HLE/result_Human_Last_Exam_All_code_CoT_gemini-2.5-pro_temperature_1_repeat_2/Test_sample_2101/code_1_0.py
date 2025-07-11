from fractions import Fraction

def solve_electron_escape_probability():
    """
    Calculates the difference in escape probability for an electron in a right isosceles triangle.
    """
    # Step 1: Explain the methodology based on a geometric probability theorem.
    print("We will solve this problem using a theorem from geometric probability.")
    print("The probability of a random ray from a random point in a triangle escaping through a given side is:")
    print("P(Side) = (angle1 + angle2) / (2 * pi), where angle1 and angle2 are the internal angles at the side's endpoints.")
    print("-" * 30)

    # Step 2: Define the angles of the isosceles right triangle in terms of fractions of pi.
    # The angles are 90 degrees (pi/2) and two 45-degree angles (pi/4).
    angle_right = Fraction(1, 2)  # Represents pi/2
    angle_acute = Fraction(1, 4)  # Represents pi/4
    print(f"The triangle's angles are pi*{angle_right}, pi*{angle_acute}, and pi*{angle_acute}.")
    print("-" * 30)

    # Step 3: Calculate the probability of escaping through the hypotenuse.
    # The hypotenuse connects the two vertices with the acute angles (pi/4).
    p_hyp_numerator = angle_acute + angle_acute
    # The probability is (numerator * pi) / (2 * pi) = numerator / 2
    p_hyp = p_hyp_numerator / 2
    
    print("1. Probability of escaping through the hypotenuse, P(hyp):")
    print(f"The angles at the endpoints are pi*{angle_acute} and pi*{angle_acute}.")
    print(f"P(hyp) = (pi*{angle_acute} + pi*{angle_acute}) / (2*pi) = (pi*{p_hyp_numerator}) / (2*pi) = {p_hyp}")
    print("-" * 30)

    # Step 4: Calculate the total probability of escaping through the two legs.
    # Each leg connects the right-angle vertex (pi/2) and an acute-angle vertex (pi/4).
    p_leg_numerator = angle_right + angle_acute
    p_leg = p_leg_numerator / 2
    p_legs_total = 2 * p_leg

    print("2. Probability of escaping through the two legs, P(legs):")
    print(f"The angles at a leg's endpoints are pi*{angle_right} and pi*{angle_acute}.")
    print(f"P(one leg) = (pi*{angle_right} + pi*{angle_acute}) / (2*pi) = (pi*{p_leg_numerator}) / (2*pi) = {p_leg}")
    print(f"Since there are two legs, P(legs) = 2 * {p_leg} = {p_legs_total}")
    print("-" * 30)

    # Step 5: Calculate the final difference.
    difference = p_hyp - p_legs_total

    print("3. The difference between the two probabilities is P(hyp) - P(legs).")
    print("\nFinal Equation:")
    print(f"{p_hyp} - {p_legs_total} = {difference}")

    # Final answer in the required format
    print(f"\n<<<{float(difference)}>>>")

solve_electron_escape_probability()