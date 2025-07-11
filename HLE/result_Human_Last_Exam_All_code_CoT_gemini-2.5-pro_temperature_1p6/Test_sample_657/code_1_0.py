import sympy as sp

def main():
    """
    This script calculates and displays the set of all proper stabilizing controllers
    for the plant H_1(s) = s / (s^2 - 1) using the Youla-Kucera parametrization.
    """
    
    # Define the symbolic variable 's' and the parameter function 'K(s)'
    s = sp.symbols('s')
    K = sp.Function('K')(s)

    # Based on the Youla-Kucera parametrization derived in the explanation,
    # we define the numerator and denominator of the controller H_2(s).
    
    # Numerator of H_2(s)
    # N_2(s) = 4(s+1)^2 + (s^2-1)K(s)
    controller_numerator = 4 * (s + 1)**2 + (s**2 - 1) * K

    # Denominator of H_2(s)
    # D_2(s) = s^2 - 1 - s*K(s)
    controller_denominator = s**2 - 1 - s * K

    print("The set of all proper stabilizing controllers H_2(s) is given in the form H_2(s) = N_2(s) / D_2(s), where K(s) is any stable and proper transfer function.")
    print("\nThe resulting parametrized transfer function is:\n")

    # Use sympy's pretty print for a clean mathematical representation
    sp.pprint(sp.Eq(sp.Function('H_2')(s), controller_numerator / controller_denominator, evaluate=False), use_unicode=True)
    
    # For a clear final answer format, we print the components explicitly
    print("\nWhere the numerator is:")
    final_numerator_expr = 4*s**2 + 8*s + 4 + (s**2 - 1)*K
    sp.pprint(final_numerator_expr)

    print("\nAnd the denominator is:")
    final_denominator_expr = s**2 - 1 - s*K
    sp.pprint(final_denominator_expr)


if __name__ == "__main__":
    main()
