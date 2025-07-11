import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron in an
    isosceles right triangle.
    """

    # --- Step 1 & 2: Define Geometry ---
    # We can set the leg length L to 1.0, as it will cancel out in the final ratio.
    leg_length = 1.0
    print(f"Let the length of the two equal legs be L = {leg_length:.1f}")

    # --- Step 3: Calculate Side Lengths ---
    hypotenuse_length = leg_length * math.sqrt(2)
    legs_total_length = leg_length + leg_length
    perimeter = legs_total_length + hypotenuse_length
    
    print(f"Length of the hypotenuse = L * sqrt(2) = {hypotenuse_length:.4f}")
    print(f"Combined length of the two legs = 2 * L = {legs_total_length:.1f}")
    print(f"Total perimeter of the triangle = 2*L + L*sqrt(2) = {perimeter:.4f}\n")

    # --- Step 4: Calculate Probabilities ---
    # P(Side) = Length(Side) / Perimeter
    prob_hypotenuse = hypotenuse_length / perimeter
    prob_legs = legs_total_length / perimeter

    print(f"Probability of escaping through the hypotenuse, P(H) = {hypotenuse_length:.4f} / {perimeter:.4f} = {prob_hypotenuse:.4f}")
    print(f"Probability of escaping through the legs, P(L) = {legs_total_length:.1f} / {perimeter:.4f} = {prob_legs:.4f}\n")

    # --- Step 5: Compute the Difference ---
    difference = prob_hypotenuse - prob_legs

    # The exact analytical solution simplifies to 2*sqrt(2) - 3
    print("The required difference is P(H) - P(L).")
    print(f"The exact expression for this difference is 2 * sqrt(2) - 3.")
    
    # As requested, printing each number in the final equation
    val_sqrt2 = math.sqrt(2)
    val_2 = 2.0
    val_3 = 3.0

    print("\nThe final equation with numerical values is:")
    print(f"{val_2:.1f} * {val_sqrt2:.4f} - {val_3:.1f} = {difference:.4f}")

    # Return the exact value as well for the final answer tag
    return 2 * math.sqrt(2) - 3

# Run the solver and print the final result.
final_answer = solve_escape_probability()
# The problem statement requests the final answer in a specific format.
# print(f"\nFinal Answer: {final_answer}")
# Example: <<<2*sqrt(2)-3>>>

if __name__ == '__main__':
    solve_escape_probability()
