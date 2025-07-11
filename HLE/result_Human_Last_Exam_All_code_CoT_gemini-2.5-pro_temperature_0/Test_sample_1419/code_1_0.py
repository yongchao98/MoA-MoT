import math

def find_fixed_point_expression():
    """
    This function calculates and prints the leading order expression for the
    Wilson-Fisher fixed point coupling u* in phi^4 theory near 4 dimensions.
    """

    # The beta function for the dimensionless coupling u in d = 4 - epsilon dimensions is:
    # beta(u) = -epsilon * u + (3 / (16 * pi^2)) * u^2

    # A fixed point u* is found by setting beta(u*) = 0. For the non-trivial
    # Wilson-Fisher fixed point, we solve:
    # epsilon = (3 / (16 * pi^2)) * u*

    # Solving for u* gives:
    # u* = (16 * pi^2 / 3) * epsilon

    # We will now print this equation, showing each number involved.
    numerator_val = 16
    denominator_val = 3
    pi_symbol = "π"
    epsilon_symbol = "ϵ"

    print("The leading order expression for the fixed point coupling u* in terms of ϵ is:")
    # The final equation is u* = C * epsilon. We specify C.
    final_equation = f"u* = ({numerator_val}{pi_symbol}² / {denominator_val}) {epsilon_symbol}"
    print(final_equation)

    # Breaking down the equation as requested:
    print("\nIn this equation:")
    print(f"  - The numerator of the coefficient involves the number '{numerator_val}' and '{pi_symbol}²'.")
    print(f"  - The denominator of the coefficient is the number '{denominator_val}'.")
    print(f"  - The expression is proportional to '{epsilon_symbol}'.")

    # Calculate the numerical value of the coefficient
    coefficient_value = (numerator_val * math.pi**2) / denominator_val

    print(f"\nThe numerical value of the coefficient ({numerator_val}π² / {denominator_val}) is approximately:")
    print(f"{coefficient_value:.4f}")

    print("\nSo, the expression can be written with the numerical approximation as:")
    print(f"u* ≈ {coefficient_value:.4f} {epsilon_symbol}")

if __name__ == "__main__":
    find_fixed_point_expression()
