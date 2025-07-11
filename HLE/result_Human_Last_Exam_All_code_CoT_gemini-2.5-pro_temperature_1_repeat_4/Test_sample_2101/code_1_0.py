import math

def solve_electron_probability():
    """
    Calculates the difference between the probability that an electron escapes through the
    hypotenuse of an isosceles right triangle and the probability that it escapes
    through either of the two legs.

    The solution is based on the principle that the escape probability for a given
    side is the ratio of that side's length to the triangle's total perimeter.
    The final simplified expression for the difference is 2*sqrt(2) - 3.
    """

    # The final simplified formula for the difference is 2*sqrt(2) - 3.
    # We will now calculate this value and show each number in the final equation.

    print("The final result is derived from the exact expression: 2 * sqrt(2) - 3")
    print("\n--- Calculation Breakdown ---")

    # First number in the equation: 2 * sqrt(2)
    sqrt_2 = math.sqrt(2)
    first_term = 2 * sqrt_2

    # Second number in the equation: -3
    second_term = -3.0

    # The final result is the sum of these two terms.
    final_difference = first_term + second_term

    print(f"The first term in the equation is 2 * sqrt(2), which is approximately: {first_term}")
    print(f"The second term in the equation is: {second_term}")
    print("\nThe final equation is the sum of these terms:")
    print(f"{first_term} + ({second_term}) = {final_difference}")

solve_electron_probability()