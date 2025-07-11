import math

def solve_electron_escape_probability():
    """
    Calculates the difference in escape probability for an electron in an isosceles right triangle.

    The solution is based on the principle from geometric probability that for a random
    point and direction in a convex shape, the escape probability through a side
    is proportional to the length of that side.
    P(escape through side S) = Length(S) / Total Perimeter
    """

    # For simplicity, we can set the length of the equal legs to 1.
    # The final result is scale-invariant.
    leg_length = 1.0

    # Calculate the length of the hypotenuse using the Pythagorean theorem
    hypotenuse_length = math.sqrt(leg_length**2 + leg_length**2)

    # Calculate the total perimeter of the triangle
    perimeter = leg_length + leg_length + hypotenuse_length

    # Calculate the probability of the electron escaping through the hypotenuse
    prob_hypotenuse = hypotenuse_length / perimeter

    # The probability of escaping through *either* leg is the sum of the probabilities
    # for each leg. Since the legs have equal length, this is 2 * P(one leg).
    prob_either_leg = (leg_length + leg_length) / perimeter

    # Calculate the final difference as requested by the problem
    difference = prob_hypotenuse - prob_either_leg
    
    # --- Outputting the results ---
    print("Step 1: Define the geometry of the isosceles right triangle.")
    print(f" - Assuming leg length = {leg_length:.1f}")
    print(f" - Calculated hypotenuse length = sqrt(2) = {hypotenuse_length:.5f}")
    print(f" - Calculated total perimeter = 2 + sqrt(2) = {perimeter:.5f}")

    print("\nStep 2: Calculate the probabilities.")
    print(f" - Probability of escaping through hypotenuse = {prob_hypotenuse:.5f}")
    print(f" - Probability of escaping through either leg = {prob_either_leg:.5f}")

    print("\nStep 3: Calculate the final difference.")
    print("   Difference = P(hypotenuse) - P(either leg)")

    # Output the final equation with the computed numbers as requested
    print("\nFinal Equation:")
    print(f"   {prob_hypotenuse:.5f} - {prob_either_leg:.5f} = {difference:.5f}")

    # The exact symbolic result is 2*sqrt(2) - 3
    exact_result = 2 * math.sqrt(2) - 3
    print(f"\nNote: The exact symbolic answer is 2*sqrt(2) - 3, which is approximately {exact_result:.5f}.")


# Run the solver
solve_electron_escape_probability()

# The final answer is the numerical result of the calculation.
final_answer = 2 * math.sqrt(2) - 3
# <<<final_answer>>>