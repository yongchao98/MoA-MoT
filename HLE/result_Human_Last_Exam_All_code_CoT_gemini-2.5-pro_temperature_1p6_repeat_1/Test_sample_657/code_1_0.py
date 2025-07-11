import sympy

def solve_controller():
    """
    This function calculates and prints the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1).
    The controller is parametrized by a function K(s).
    """

    # Define the symbolic variable 's' and the parameter function K(s)
    s = sympy.Symbol('s')
    K_func = sympy.Function('K')(s)

    # The components of the Youla-Kucera parametrization derived from the plan
    # Numerator_H2 = 4(s+1)^2 + (s^2-1)K(s)
    # Denominator_H2 = (s^2-1) - sK(s)

    # Define the numerator and denominator of the controller H_2(s)
    num_H2 = 4 * (s + 1)**2 + (s**2 - 1) * K_func
    den_H2 = (s**2 - 1) - s * K_func

    # Expand the expressions for clarity
    expanded_num = sympy.expand(num_H2)
    expanded_den = sympy.expand(den_H2)

    # Print the results
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("H_2(s) = N(s) / D(s)\n")

    print("Numerator N(s):")
    sympy.pprint(expanded_num)
    print("\nDenominator D(s):")
    sympy.pprint(expanded_den)

    print("\nNote: For the controller H_2(s) to be proper and stabilizing, the parameter K(s) must be chosen")
    print("as any stable and proper rational transfer function.")

# Execute the function to find and print the solution
solve_controller()