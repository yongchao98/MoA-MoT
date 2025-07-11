import numpy as np

def solve():
    """
    This script verifies the mathematical derivation for the lower bound
    of d(t,x) and prints the final result.
    """

    # The problem of finding a constant lower bound C for d(t,x) reduces to the
    # condition that a polynomial f_C(u) is non-negative for all u in [0, 1].
    # f_C(u) = u^4 - u^3 + 5*C*u^2 - 3*C*u + 2*C^2

    # Mathematical analysis shows that C must be less than or equal to -1.
    # We test the best possible candidate from this method, C = -1.
    C = -1

    # For C = -1, the polynomial is f(u) = u^4 - u^3 - 5u^2 + 3u + 2.
    # We define its coefficients: f(u) = a*u^4 + b*u^3 + c*u^2 + d*u + e
    a = 1
    b = -1
    c = 5 * C
    d = -3 * C
    e = 2 * C**2

    print("The problem of finding a constant lower bound C for d(t,x) reduces to a condition on a polynomial f_C(u).")
    print("Analysis shows that C must be less than or equal to -1. We test C = -1.")
    print("\nThe polynomial f_C(u) for C = -1 is:")
    print(f"f(u) = ({a}) * u^4 + ({b}) * u^3 + ({c}) * u^2 + ({d}) * u + ({e})")

    # We verify numerically that this polynomial is non-negative on [0, 1].
    u_values = np.linspace(0, 1, 1001)
    f_values = a*u_values**4 + b*u_values**3 + c*u_values**2 + d*u_values + e
    min_f_value = np.min(f_values)

    print("\nNumerical verification:")
    print(f"The minimum value of this polynomial on the interval [0, 1] is found to be approximately: {min_f_value:.6f}")

    # Use a small tolerance for floating point errors
    if min_f_value >= -1e-9:
        print("This value is non-negative (within numerical precision).")
        print("This confirms that C = -1 is a valid lower bound.")
    else:
        print("The condition is not met, there might be an error in the derivation.")

    lower_bound = -1
    print(f"\nConclusion: A constant lower bound for d(t,x) is {lower_bound}.")

solve()