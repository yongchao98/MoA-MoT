import sympy

def solve_and_print_expectation():
    """
    Calculates the expected value of T symbolically using sympy and prints the result.
    The expected value is E[T] = sum_{i=1 to infinity} i / (3^i * (i+1)!).
    """
    # Define the symbolic variable for the summation index
    i = sympy.Symbol('i', integer=True, positive=True)

    # Define the term in the summation
    term = i / (3**i * sympy.factorial(i + 1))

    # Calculate the infinite sum using sympy
    expected_value_expr = sympy.summation(term, (i, 1, sympy.oo))

    # The result is of the form a - b * exp(c/d). We extract these numbers.
    # The expression is 3 - 2*exp(1/3)
    # expr.args gives (3, -2*exp(1/3))
    a = expected_value_expr.args[0]
    
    # term2 is -2*exp(1/3). term2.args is (-2, exp(1/3))
    term2 = expected_value_expr.args[1]
    b_neg = term2.args[0]
    b = -b_neg
    
    # exp_term is exp(1/3). exp_term.args is (1/3,)
    exp_term = term2.args[1]
    power_frac = exp_term.args[0] # This is a sympy Fraction
    c = power_frac.p
    d = power_frac.q

    # Print the final equation and its components
    print(f"The expected value of T is given by the equation: {a} - {b} * e^({c}/{d})")
    print("The numbers in the final equation are:")
    print(a)
    print(b)
    print(c)
    print(d)

    # For verification, print the symbolic expression from sympy and its numerical value
    print("\nSymbolic result from sympy:")
    sympy.pprint(expected_value_expr, use_unicode=True)
    
    print(f"\nNumerical value: {expected_value_expr.evalf()}")

solve_and_print_expectation()