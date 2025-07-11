import math

def solve():
    """
    This function solves the problem based on the analytical derivation.
    """
    # Given parameters from the problem description
    # |V| = n
    n = 99

    # The problem asks to compute the sum S = sum_{w in V^n} a(w),
    # where a(w) = (n + 1 - |unique_tokens(w)|)^(-1).
    # The analytical derivation of this sum leads to a simple closed-form solution:
    # S = (n + 1)^(n - 1)

    # Let's plug in the value of n=99 to calculate the sum.
    base = n + 1
    exponent = n - 1

    # The result is base^exponent
    # For n=99, this is 100^98.

    # The question asks to write the answer as a power of 10.
    # We can rewrite the base as a power of 10: 100 = 10^2.
    # So, Sum = (10^2)^98 = 10^(2 * 98) = 10^196.
    power_of_10_exponent = 2 * exponent

    # Print the steps of the final calculation as requested.
    print(f"The vocabulary size is n = {n}.")
    print("The formula for the sum is (n + 1)^(n - 1).")
    print("\nSubstituting the value of n:")
    print(f"Sum = ({n} + 1)^({n} - 1)")
    # The output string explicitly shows each number in the equation.
    print(f"Sum = {base}^{exponent}")
    
    print("\nExpressing the result as a power of 10:")
    print(f"Sum = (10^2)^{exponent} = 10^(2 * {exponent}) = 10^{power_of_10_exponent}")

solve()