import decimal

def solve_for_T():
    """
    Calculates the time T based on the solvability condition of the boundary-value problem.
    """
    # Set the precision for decimal calculations. The exponents involved require
    # a high precision to maintain accuracy. The largest exponent is in the order
    # of 5,010,000, so we set the precision to a slightly higher value.
    decimal.getcontext().prec = 5010010

    # Define the given parameters as Decimal objects to handle the extreme magnitudes.
    alpha = decimal.Decimal('1e10000')
    x0 = decimal.Decimal('1e-5000000')
    
    # y0 is given to be equal to x0.
    y0 = x0

    # The solvability condition is CT = alpha, where C is derived from the
    # first-order perturbation terms.
    # C = x0/(1-x0) + y0/(1-y0)
    # Since x0 = y0, C = 2 * x0 / (1 - x0)
    #
    # We solve for T: T = alpha / C
    # T = alpha * (1 - x0) / (2 * x0)

    one = decimal.Decimal(1)
    two = decimal.Decimal(2)

    # Calculate T using the derived formula.
    T = alpha * (one - x0) / (two * x0)

    # The problem asks to output each number in the final equation.
    # The final equation is T * (2 * x0 / (1 - x0)) = alpha.
    # We print the main variables and the final result for T.
    
    print("Derived final equation: T * (2 * x0 / (1 - x0)) = alpha")
    print("-" * 30)
    print(f"Value of alpha: {alpha}")
    print(f"Value of x0: {x0}")
    print("-" * 30)
    print(f"Calculated value of T:")
    # The string representation of T is very long, but it is the exact result
    # of the high-precision calculation.
    print(T)

solve_for_T()