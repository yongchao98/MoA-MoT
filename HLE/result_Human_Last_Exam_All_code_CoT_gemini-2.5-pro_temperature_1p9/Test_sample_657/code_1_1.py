import sympy

def main():
    """
    This function prints the set of all proper stabilizing controllers H2(s)
    for the given system H1(s).
    """

    # Define the symbolic variable s and the controller H2(s) for printing
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Construct the numerator and denominator of the controller transfer function
    # Numerator: -(s+1) * (2 + (s-1)*K(s))
    num_part_1_coeff = 2
    num_part_2_factor_coeff = 1
    num_part_2_factor = s - num_part_2_factor_coeff
    num_inner = num_part_1_coeff + num_part_2_factor * K
    numerator = -(s + 1) * num_inner

    # Denominator: 1 + s*K(s)
    den_part_1_coeff = 1
    denominator = den_part_1_coeff + s * K

    # Pretty print the final expression for the controller H_2(s)
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    h2_expression = sympy.Eq(sympy.Function('H_2')(s), numerator / denominator)
    sympy.pprint(h2_expression, use_unicode=True)

    # Print the conditions on the parameter K(s)
    print("\nwhere K(s) is any stable rational function that satisfies both of the following conditions:")
    print("1. K(s) must be strictly proper, i.e., lim K(s) = 0 as s -> infinity.")
    print("2. The limit of s*K(s) as s -> infinity must be -2.")
    print("\nAn example of such a function is K(s) = -2 / (s + a) for any a > 0.")

if __name__ == '__main__':
    main()
