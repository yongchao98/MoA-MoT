import math

def get_approximation_formula():
    """
    Calculates and prints the parameters for the analytical approximation of the integral
    for the small epsilon regime. The formula is of the form I(epsilon) = C * epsilon^k.
    """

    # The approximation is based on the lowest power term in the denominator's polynomial,
    # which is 9.0 * x^5.0. Thus, n=5.0 and the coefficient is c_n=9.0.
    n = 5.0
    c_n = 9.0

    # The exponent k in the approximation I(epsilon) ~ C * epsilon^k is given by k = (1 - n) / n.
    k = (1.0 - n) / n

    # The coefficient C is given by the formula C = pi / (n * c_n^(1/n) * sin(pi/n)).
    C = math.pi / (n * math.pow(c_n, 1.0 / n) * math.sin(math.pi / n))

    print("The analytical formula for the integral I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) = C * epsilon^k")
    print("\nThe values for the numbers in this equation are:")
    print(f"C = {C}")
    print(f"k = {k}")

    print("\nTherefore, the final approximate formula is:")
    print(f"I(epsilon) = {C:.5f} * epsilon^({k:.2f})")

get_approximation_formula()