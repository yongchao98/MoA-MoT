import math

def solve():
    """
    Calculates the constants for the analytical approximation of the integral I(epsilon).
    The approximation has the form I(epsilon) ≈ C * epsilon^(-p).
    """
    # Parameters from the dominant term in the denominator polynomial, a*x^p
    a = 9.0
    n = 5.0

    # The exponent of epsilon in the approximation
    p = (n - 1) / n

    # The coefficient C in the approximation
    # C = (pi / n) / (a^(1/n) * sin(pi/n))
    C = (math.pi / n) / (a**(1/n) * math.sin(math.pi / n))

    print("The analytical formula for I(epsilon) in the small epsilon regime is of the form:")
    print("I(epsilon) ≈ C * epsilon^(-p)")
    print("\nCalculated values for the constants:")
    print(f"The coefficient C is: {C}")
    print(f"The exponent p is: {p}")

    print("\nThus, the final approximate formula is:")
    print(f"I(epsilon) ≈ {C} * epsilon^(-{p})")

solve()