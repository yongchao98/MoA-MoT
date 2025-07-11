import math

def solve_electron_escape_probability():
    """
    Calculates the difference between the probability that an electron escapes
    through the hypotenuse and the probability that it escapes through either
    of the two legs of an isosceles right triangle.

    The solution is based on the principle that the probability of escape
    through a given side is proportional to the length of that side compared
    to the total perimeter.
    """

    # Let the length of the two equal legs be L. For simplicity, we can set L=1,
    # as it will cancel out in the probability ratios.
    L = 1.0
    
    # 1. Calculate the lengths of the sides
    leg1_length = L
    leg2_length = L
    hypotenuse_length = math.sqrt(leg1_length**2 + leg2_length**2)

    # 2. Calculate the total perimeter
    perimeter = leg1_length + leg2_length + hypotenuse_length

    # 3. Calculate the probability of escaping through each side
    prob_hypotenuse = hypotenuse_length / perimeter
    prob_leg1 = leg1_length / perimeter
    prob_leg2 = leg2_length / perimeter

    # 4. Calculate the required difference
    difference = prob_hypotenuse - (prob_leg1 + prob_leg2)

    # 5. Print the step-by-step derivation of the symbolic answer
    print("Step-by-step derivation:")
    print("Let the length of the legs be L. We can set L=1 for simplicity.")
    print(f"Length of leg 1 = {leg1_length}")
    print(f"Length of leg 2 = {leg2_length}")
    print(f"Length of hypotenuse = sqrt({leg1_length}^2 + {leg2_length}^2) = sqrt(2)")
    print(f"Perimeter = {leg1_length} + {leg2_length} + sqrt(2) = 2 + sqrt(2)")
    print("-" * 30)
    print("The probability of escape through a side is its length divided by the perimeter.")
    print("P(hypotenuse) = sqrt(2) / (2 + sqrt(2))")
    print("P(leg) = 1 / (2 + sqrt(2))")
    print("-" * 30)
    print("The required difference is P(hypotenuse) - (P(leg1) + P(leg2))")
    print("= [sqrt(2) / (2 + sqrt(2))] - [1 / (2 + sqrt(2)) + 1 / (2 + sqrt(2))]")
    print("= [sqrt(2) / (2 + sqrt(2))] - [2 / (2 + sqrt(2))]")
    print("= (sqrt(2) - 2) / (2 + sqrt(2))")
    print("\nTo simplify, multiply numerator and denominator by the conjugate (2 - sqrt(2)):")
    print("= [(sqrt(2) - 2) * (2 - sqrt(2))] / [(2 + sqrt(2)) * (2 - sqrt(2))]")
    print("= [2*sqrt(2) - 2 - 4 + 2*sqrt(2)] / [4 - 2]")
    print("= (4*sqrt(2) - 6) / 2")
    print("= 2*sqrt(2) - 3")
    print("-" * 30)
    print(f"Numerical value: {difference}")

solve_electron_escape_probability()
