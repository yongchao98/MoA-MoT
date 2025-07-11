import math
from fractions import Fraction

def solve_electron_escape_probability():
    """
    Solves for the difference in probabilities using the integral geometry theorem.
    <α_s1 + α_s2> = π + β
    """

    # The problem asks for P(hyp) - P(legs), where P(legs) = P(leg1) + P(leg2).
    # This is equivalent to (<α_hyp> - <α_leg1> - <α_leg2>) / (2π).

    # The internal angles of the isosceles right triangle are:
    # beta_right_angle = π/2
    # beta_acute_angle = π/4

    # From the theorem, we can set up a system of equations for the average
    # subtended angles. Let the variables represent the average angles divided by π.
    # x = <α_leg1> / π
    # y = <α_leg2> / π
    # z = <α_hyp> / π
    #
    # 1) Legs meet at the right angle (π/2): x + y = 1 + 1/2 = 3/2
    # 2) Leg1 and Hypotenuse meet at an acute angle (π/4): x + z = 1 + 1/4 = 5/4
    # 3) Leg2 and Hypotenuse meet at an acute angle (π/4): y + z = 1 + 1/4 = 5/4

    # From (2) and (3), we can see that x must be equal to y.
    # Substitute x=y into (1):
    # 2x = 3/2  =>  x = 3/4
    # So, avg_alpha_leg1 / π = 3/4 and avg_alpha_leg2 / π = 3/4.
    avg_alpha_leg_coeff = Fraction(3, 4)

    # Substitute x = 3/4 into (2) to find z:
    # 3/4 + z = 5/4  => z = 5/4 - 3/4 = 2/4 = 1/2
    avg_alpha_hyp_coeff = Fraction(1, 2)

    # Now calculate the probabilities P(side) = <α_side> / (2π) = (<α_side>/π) / 2
    prob_leg1 = avg_alpha_leg_coeff / 2
    prob_leg2 = avg_alpha_leg_coeff / 2
    prob_hyp = avg_alpha_hyp_coeff / 2

    # Total probability for escaping through either of the two legs
    prob_legs = prob_leg1 + prob_leg2

    # The difference between the probability of escaping through the hypotenuse
    # and the probability of escaping through the legs.
    difference = prob_hyp - prob_legs

    print("--- Probabilities ---")
    print(f"Probability of escaping through the hypotenuse: {prob_hyp.numerator}/{prob_hyp.denominator}")
    print(f"Probability of escaping through the legs: {prob_legs.numerator}/{prob_legs.denominator}")
    print("\n--- Final Calculation ---")
    print(f"The difference is the probability of escaping through the hypotenuse minus the probability of escaping through the legs.")
    print(f"Difference = P(hypotenuse) - P(legs)")
    print(f"Difference = {prob_hyp.numerator}/{prob_hyp.denominator} - {prob_legs.numerator}/{prob_legs.denominator}")
    print(f"Difference = {difference.numerator}/{difference.denominator}")
    print(f"Difference = {float(difference)}")

solve_electron_escape_probability()
<<< -0.5 >>>