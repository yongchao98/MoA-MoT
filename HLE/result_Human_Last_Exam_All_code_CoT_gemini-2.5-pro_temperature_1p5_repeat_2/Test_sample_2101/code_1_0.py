import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron
    in an isosceles right triangular material.
    """

    # Let the length of the two equal legs be L. We can set L=1 for simplicity.
    L = 1.0

    # Calculate the lengths of the sides
    length_leg = L
    length_hypotenuse = L * math.sqrt(2)

    # Calculate the total perimeter of the triangle
    perimeter = 2 * length_leg + length_hypotenuse

    # The probability of escaping through a side is proportional to its length.
    # P(side) = length_of_side / perimeter
    prob_hypotenuse = length_hypotenuse / perimeter
    prob_one_leg = length_leg / perimeter

    # The total probability of escaping through the legs is the sum of the probabilities for each leg.
    prob_legs_total = 2 * prob_one_leg

    # Calculate the required difference
    difference = prob_hypotenuse - prob_legs_total

    # --- Output the results step-by-step ---
    print("Let the length of the legs be L.")
    print("The length of the hypotenuse is L * sqrt(2).")
    print(f"Thus, the perimeter is 2*L + L*sqrt(2) = L * (2 + sqrt(2)).")
    print("-" * 30)

    # Print the equation for the probability of escaping through the hypotenuse
    print("Probability of escaping through hypotenuse:")
    print(f"P(hypotenuse) = (L * sqrt(2)) / (L * (2 + sqrt(2))) = sqrt(2) / (2 + sqrt(2))")
    print(f"P(hypotenuse) ~= {prob_hypotenuse:.4f}")
    print("-" * 30)
    
    # Print the equation for the probability of escaping through the legs
    print("Total probability of escaping through the two legs:")
    print(f"P(legs) = (2 * L) / (L * (2 + sqrt(2))) = 2 / (2 + sqrt(2))")
    print(f"P(legs) ~= {prob_legs_total:.4f}")
    print("-" * 30)

    # Print the final calculation for the difference
    print("The difference is P(hypotenuse) - P(legs):")
    print(f"Difference = (sqrt(2) / (2 + sqrt(2))) - (2 / (2 + sqrt(2)))")
    print(f"Difference = (sqrt(2) - 2) / (2 + sqrt(2))")
    print("-" * 30)

    print("Simplifying the expression gives 2*sqrt(2) - 3.")
    print(f"Final calculated difference: {difference}")

if __name__ == "__main__":
    solve_escape_probability()
    final_answer = 2 * math.sqrt(2) - 3
    # The final answer is wrapped according to the format specification
    print(f'<<<{final_answer}>>>')
