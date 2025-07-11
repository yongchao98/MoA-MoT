import math

def calculate_max_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for the described rotating charged sphere.

    The derivation shows the final ratio is given by the expression:
    Ratio = (pi / 3) * (16 - 7 * sqrt(2))
    """

    # Define the constants from the derived formula
    pi = math.pi
    sqrt2 = math.sqrt(2)

    # The terms in the equation
    term_16 = 16
    term_7_sqrt_2 = 7 * sqrt2
    factor_pi_over_3 = pi / 3

    # Calculate the final result
    max_ratio = factor_pi_over_3 * (term_16 - term_7_sqrt_2)

    # Print the explanation and the result
    print("The problem is to find the maximum ratio of bidirectional conical power to line intensity.")
    print("The derivation leads to the following analytical expression for the ratio:")
    print("Ratio = (π / 3) * (16 - 7 * √2)\n")

    print("Let's calculate this value step-by-step:")
    print(f"π (pi) ≈ {pi}")
    print(f"√2 (sqrt(2)) ≈ {sqrt2}")
    print(f"The term (16 - 7 * √2) evaluates to: {term_16} - {term_7_sqrt_2} = {term_16 - term_7_sqrt_2}")
    print(f"The factor (π / 3) evaluates to: {pi} / 3 = {factor_pi_over_3}\n")

    print("Final Result:")
    print(f"Ratio = {factor_pi_over_3} * ({term_16 - term_7_sqrt_2})")
    print(f"Ratio ≈ {max_ratio}")

if __name__ == "__main__":
    calculate_max_ratio()