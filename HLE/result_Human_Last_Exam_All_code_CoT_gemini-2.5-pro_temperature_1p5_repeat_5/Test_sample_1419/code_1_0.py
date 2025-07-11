import math

def display_fixed_point_expression():
    """
    This function prints the leading order expression for the fixed point
    coupling u* in phi^4 theory near four dimensions.
    """

    # The expression for the fixed point coupling u* is (16 * pi^2 / 3) * epsilon.
    # We define the numerical coefficients here as requested.
    numerator_coefficient = 16
    denominator_coefficient = 3
    pi_symbol = "pi"
    epsilon_symbol = "epsilon"

    print("The leading order expression for the fixed point coupling u* in phi^4 theory is given by the Wilson-Fisher fixed point.")
    print("The equation is:")
    # Using 'pi^2' for clarity, as math.pi**2 would be a numerical evaluation.
    print(f"u* = ({numerator_coefficient} * {pi_symbol}^2 / {denominator_coefficient}) * {epsilon_symbol}")

if __name__ == "__main__":
    display_fixed_point_expression()