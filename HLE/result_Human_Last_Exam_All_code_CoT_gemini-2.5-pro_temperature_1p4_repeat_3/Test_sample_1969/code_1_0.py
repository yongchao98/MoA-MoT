def print_shapley_value_formula():
    """
    This function prints the derived formula for the Shapley value c_k.
    The problem asks for the formula for c_k for a fixed n.
    The formula is: c_k = k * S1 * (S1^2 - k*S1 + S2)
    where S1 and S2 are sums of powers of the first n integers.
    """

    # The problem asks to output each number in the final equation.
    # The coefficients derived for the final simplified formula are all 1.
    k_coeff = 1
    S1_coeff1 = 1
    S1_sq_coeff = 1
    k_S1_coeff = 1
    S2_coeff = 1

    formula_str = (
        f"The exact amount of money c_k that person p_k gets is given by the formula:\n\n"
        f"c_k = {k_coeff} * k * S1 * ({S1_sq_coeff} * S1^2 - {k_S1_coeff} * k * S1 + {S2_coeff} * S2)\n\n"
        f"where:\n"
        f"n is the total number of people,\n"
        f"k is the index of the person (from 1 to n),\n"
        f"S1 is the sum of the first n integers: S1 = 1 + 2 + ... + n = n*(n+1)/2,\n"
        f"S2 is the sum of the first n squared integers: S2 = 1^2 + 2^2 + ... + n^2 = n*(n+1)*(2*n+1)/6."
    )

    print(formula_str)

print_shapley_value_formula()