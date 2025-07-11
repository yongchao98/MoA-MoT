import math

def calculate_expected_max_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the provided distributions.

    The problem involves:
    1. Yearly maximum magnitude X ~ Pareto(alpha=2, xm=1)
       - CDF: F_X(x) = 1 - 1/x^2
    2. Number of years N ~ LogSeries(p=1/2)
       - PGF: G_N(z) = log(1 - p*z) / log(1 - p)

    The CDF of the overall maximum M = max(X_1, ..., X_N) is F_M(x) = G_N(F_X(x)).
    After substitution and simplification:
    F_M(x) = (log(2*x^2) - log(x^2 + 1)) / log(2)

    The expectation E[M] is the integral of the survival function from 1 to infinity:
    E[M] = Integral from 1 to inf of [1 - F_M(x)] dx
    1 - F_M(x) = log(1 + 1/x^2) / log(2)

    The integral of log(1 + 1/x^2) from 1 to infinity evaluates to (pi/2 - log(2)).
    So, E[M] = (pi/2 - log(2)) / log(2)
             = pi / (2 * log(2)) - 1
    """

    # The final equation is: E = (pi/2 - log(2)) / log(2)
    # We will print the values of the components of this equation.

    # Value of pi
    val_pi = math.pi
    
    # Value of the natural logarithm of 2
    val_log_2 = math.log(2)

    # Calculate the numerator: pi/2 - log(2)
    numerator = (val_pi / 2) - val_log_2

    # The denominator is just log(2)
    denominator = val_log_2

    # Calculate the final expected value
    expected_value = numerator / denominator

    print("This script calculates the expected maximum earthquake magnitude.")
    print("The analytical formula for the expectation is: (pi/2 - log(2)) / log(2)\n")
    print("--- Components of the Final Equation ---")
    print(f"Value of pi (π): {val_pi}")
    print(f"Value of natural log of 2 (log(2)): {val_log_2}")
    print(f"Numerator (pi/2 - log(2)): {numerator}")
    print(f"Denominator (log(2)): {denominator}\n")
    print("--- Final Result ---")
    print(f"The expected maximum earthquake magnitude is: {expected_value}")

calculate_expected_max_magnitude()
<<<1.2661801235332219>>>