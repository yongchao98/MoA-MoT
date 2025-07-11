import math

def solve_limit():
    """
    Calculates the limit of n * P(n) as n -> infinity.

    The derivation shows the limit converges to the value (2 * sqrt(3)) / pi.
    This script calculates and prints this value, showing the components of the final fraction.
    """

    # Numerator of the final expression is 2 * sqrt(3)
    numerator_coeff = 2
    sqrt_base = 3

    # Denominator is pi
    denominator_val = math.pi

    # Calculate the numerator value
    numerator_val = numerator_coeff * math.sqrt(sqrt_base)

    # Calculate the final result of the limit
    limit_value = numerator_val / denominator_val

    # Print the explanation and the result as per the derivation
    print("The limit of n*P(n) is given by the expression: (2 * sqrt(3)) / pi")
    print("\n--- Equation Components ---")
    print(f"Numerator coefficient: {numerator_coeff}")
    print(f"Number under square root: {sqrt_base}")
    print(f"Denominator (pi): {denominator_val}")
    print("\n--- Final Equation ---")
    print(f"({numerator_coeff} * sqrt({sqrt_base})) / {denominator_val}")
    print("\n--- Result ---")
    print(f"The calculated limit is: {limit_value}")

solve_limit()