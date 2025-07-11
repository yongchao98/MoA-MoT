import math

def display_lower_bound():
    """
    This function prints the mathematical formula for the lower bound on E[S]
    as requested by the user.
    """

    # Symbolic representations for the variables in the formula
    n_symbol = 'n'
    alpha_symbol = 'alpha'
    pi_symbol = 'pi'

    # The constants in the derived expression
    constant_value = 2
    divisor = 6
    
    # Constructing the formula string. The problem requires outputting each number
    # in the final equation.
    bound_formula = (
        f"{n_symbol} * {alpha_symbol} + "
        f"{n_symbol} * ({constant_value} - {pi_symbol}^2 / {divisor})"
    )

    print("A lower bound on E[S] is given by the expression:")
    print(bound_formula)

if __name__ == '__main__':
    display_lower_bound()