import sympy

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth moment of conductance to its average
    for a disordered Majorana wire at the critical point.
    """
    # Define the symbolic variables
    # n: moment order (positive integer)
    # X: upper bound of the uniform distribution for y = -ln(g), proportional to wire length L
    # y: the variable of integration, y = -ln(g)
    n = sympy.Symbol('n', positive=True, integer=True)
    X = sympy.Symbol('X', positive=True)
    y = sympy.Symbol('y', positive=True)

    # The probability distribution of y is P(y) = 1/X for y in [0, X].
    # The integrand for calculating the n-th moment of g = exp(-y) is g^n * P(y).
    integrand = sympy.exp(-n * y) / X

    # Step 3: Calculate the exact n-th moment <g^n> as a function of X
    moment_n_exact = sympy.integrate(integrand, (y, 0, X))
    # This evaluates to (1 - exp(-n*X)) / (n*X)

    # We need the ratio <g^4> / <g>
    moment_4_exact = moment_n_exact.subs(n, 4)
    moment_1_exact = moment_n_exact.subs(n, 1)

    ratio_exact = moment_4_exact / moment_1_exact

    # Step 4: Evaluate the ratio in the limit of large wire length (X -> infinity)
    final_ratio = sympy.limit(ratio_exact, X, sympy.oo)

    # Step 5: Print the numbers in the final equation as requested.
    # The final equation is Ratio = Numerator / Denominator
    if final_ratio.is_Rational:
        numerator = final_ratio.p
        denominator = final_ratio.q
        
        print("The ratio between the fourth statistical moment of the dimensionless conductance and its average value is calculated.")
        print("Final Equation: Ratio = <g^4> / <g>")
        print(f"The numerator of the final ratio is: {numerator}")
        print(f"The denominator of the final ratio is: {denominator}")
        print(f"The value of the ratio is: {numerator}/{denominator} = {float(final_ratio)}")
    else:
        print(f"The value of the ratio is: {final_ratio}")

    return final_ratio

if __name__ == "__main__":
    ratio = solve_conductance_ratio()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<{ratio}>>>")
