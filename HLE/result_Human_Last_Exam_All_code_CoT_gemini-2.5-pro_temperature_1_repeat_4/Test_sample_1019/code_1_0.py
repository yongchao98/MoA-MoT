import math

def print_hamiltonicity_threshold_formula():
    """
    This function prints the formula for the d-threshold for Hamiltonicity
    based on the provided parameters. The formula is symbolic as n and eta
    are not given specific values.
    """

    # The only explicit integer constant in the numerator of the formula is 1.
    # As requested, we will output this number within the final equation.
    constant_one = 1

    # The formula for the threshold p is constructed as a string.
    # We use 'omega(1)' to represent a function that tends to infinity with n.
    # This is standard notation in random graph theory.
    formula = (
        f"p = (log(n) + (eta - {constant_one}) * log(log(n)) + omega({constant_one})) / n"
    )

    print("The d-threshold for Hamiltonicity in the given range is:")
    print(formula)

print_hamiltonicity_threshold_formula()
