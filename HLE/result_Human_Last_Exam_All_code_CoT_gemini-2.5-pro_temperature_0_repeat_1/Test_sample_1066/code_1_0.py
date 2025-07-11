import math

def print_lower_bound_expression():
    """
    This function prints the formula for the lower bound on the expected detection score E[S].
    The formula involves the number of tokens 'n', the average entropy 'alpha',
    and the mathematical constant 'pi'.
    """

    # The variables 'n' and 'alpha' are given symbolically. We represent them as strings.
    n_symbol = "n"
    alpha_symbol = "alpha"

    # The constant 'pi' is a number in the equation.
    pi_value = math.pi

    # The final expression for the lower bound is n * (alpha - ln(pi)).
    # We print the equation, showing each symbolic part and the constant pi.
    print("The lower bound on E[S] is given by the expression:")
    print(f"{n_symbol} * ({alpha_symbol} - ln({pi_value}))")

if __name__ == "__main__":
    print_lower_bound_expression()