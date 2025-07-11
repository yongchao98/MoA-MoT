import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron in an
    isosceles right triangle.
    """
    # The physical principle is that for a uniform and isotropic source inside a
    # convex shape, the escape probability through a side is proportional to its length.
    # P(side) = Length(side) / Total_Perimeter

    # Let L be the length of the two equal legs. We can use L=1 as it cancels out.
    L = 1.0
    print("Step 1: Define the geometry of the triangle.")
    print(f"Let the length of the two equal legs be L = {L:.1f}")

    # Calculate side lengths and perimeter
    length_leg = L
    length_hypotenuse = L * math.sqrt(2)
    perimeter = 2 * length_leg + length_hypotenuse

    print(f"The length of each leg is {length_leg:.4f}")
    print(f"The length of the hypotenuse is sqrt({L:.1f}^2 + {L:.1f}^2) = {length_hypotenuse:.4f}")
    print(f"The total perimeter is {length_leg:.4f} + {length_leg:.4f} + {length_hypotenuse:.4f} = {perimeter:.4f}\n")

    # Calculate probabilities
    prob_hypotenuse = length_hypotenuse / perimeter
    prob_leg = length_leg / perimeter
    prob_both_legs = 2 * prob_leg

    print("Step 2: Calculate the probabilities based on side lengths.")
    print(f"Probability of escape through the hypotenuse (P_h) = Length_hypotenuse / Perimeter")
    print(f"P_h = {length_hypotenuse:.4f} / {perimeter:.4f} = {prob_hypotenuse:.4f}\n")

    print(f"Probability of escape through one leg (P_l) = Length_leg / Perimeter")
    print(f"P_l = {length_leg:.4f} / {perimeter:.4f} = {prob_leg:.4f}\n")

    print("Step 3: Calculate the requested difference.")
    print("The difference is between the probability of escaping through the hypotenuse")
    print("and the probability of escaping through either of the two legs (which is 2 * P_l).\n")
    print("Difference = P_h - 2 * P_l")

    # Calculate the final numerical result
    final_difference = prob_hypotenuse - prob_both_legs

    print(f"Difference = {prob_hypotenuse:.4f} - 2 * {prob_leg:.4f}")
    print(f"Difference = {prob_hypotenuse:.4f} - {prob_both_legs:.4f} = {final_difference:.4f}\n")

    # The exact symbolic result is 2*sqrt(2) - 3.
    # We will print the final equation with its components.
    num1 = 2
    num2 = 2
    num3 = 3
    exact_result = num1 * math.sqrt(num2) - num3

    print("Step 4: State the exact result from symbolic simplification.")
    print("The simplified exact expression for the difference is 2*sqrt(2) - 3.")
    print("Final equation and result:")
    print(f"{num1} * sqrt({num2}) - {num3} = {exact_result}")

solve_escape_probability()