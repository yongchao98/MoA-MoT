import math

def solve():
    """
    Calculates the constants for the analytical approximation of the integral I(epsilon).
    The approximation has the form I(epsilon) approx C * epsilon^p.
    """

    # From the analysis, the dominant term in the denominator for small x is 9.0 * x^5.0.
    # We denote this as c * x^k.
    c = 9.0
    k = 5.0

    # The exponent p is given by the formula (1-k)/k.
    p = (1.0 - k) / k

    # The coefficient C is given by the formula c^(-1/k) * integral(1/(1+u^k)du).
    # The value of the integral is (pi/k) / sin(pi/k).
    integral_part = (math.pi / k) / math.sin(math.pi / k)
    C = math.pow(c, -1.0/k) * integral_part

    # Print the derived formula and the calculated numerical values.
    print("The analytical formula that approximates I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) ≈ C * epsilon^p")
    print(f"\nBased on the leading term {c}*x^{k}, the constants are derived as:")
    print(f"p = (1 - {k}) / {k}")
    print(f"C = {c}^(-1/{k}) * (pi/{k}) / sin(pi/{k})")
    
    print("\nNumerically, these values are:")
    print(f"p = {p}")
    print(f"C = {C}")
    
    print("\nThus, the final approximate formula is:")
    # The final equation with all numbers included, as requested.
    print(f"I(epsilon) ≈ {C} * epsilon^({p})")

solve()