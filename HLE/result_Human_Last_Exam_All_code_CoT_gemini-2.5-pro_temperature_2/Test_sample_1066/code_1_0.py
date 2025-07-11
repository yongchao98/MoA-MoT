import math

def display_lower_bound_formula():
    """
    This function prints the derived lower bound for the expected detection statistic E[S].
    The derivation mentioned in the plan leads to a formula involving the number of tokens (n),
    the average entropy per token (alpha), and the constants pi and ln(2).
    """

    # Symbolic variables in the formula
    n = "n"
    alpha = "alpha"

    # Numerical values of the constants in the formula
    pi_value = math.pi
    pi_squared_over_6 = pi_value**2 / 6
    ln_of_2 = math.log(2)

    print("Based on the derivation, a lower bound on the expected detection statistic E[S] is:")
    
    # Using symbolic names for constants for clarity
    bound_formula_symbolic = f"{n} * (pi**2 / 6 - {alpha} / ln(2))"
    print(bound_formula_symbolic)
    
    print("\nWhere the numerical values for the constants are:")
    print(f"pi = {pi_value}")
    print(f"pi**2 / 6 = {pi_squared_over_6}")
    print(f"ln(2) = {ln_of_2}")

    # Showing the formula with the calculated constant values for easier use.
    # Note: 'n' and 'alpha' remain as variables to be provided for a specific text.
    print(f"\nFinal numerical form of the equation for E[S]'s lower bound:")
    final_bound = f"{n} * ({pi_squared_over_6} - {alpha} / {ln_of_2})"
    print(final_bound)


if __name__ == '__main__':
    display_lower_bound_formula()
