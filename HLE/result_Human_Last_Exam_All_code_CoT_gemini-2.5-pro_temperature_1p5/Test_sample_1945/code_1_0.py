import math

def calculate_probability(n):
    """
    Calculates the probability Pr(n) for a given n.

    Pr(n) is the probability that a randomly selected coprime pair (p,q)
    with 1 < p,q <= n is a "good pair".
    A pair is "good" if p + q <= n + 1.
    """
    print(f"Calculating for n = {n}")
    
    # We are considering pairs (p,q) such that 1 < p <= n and 1 < q <= n.
    # The condition for a good pair is p + q <= n + 1.

    total_coprime_pairs = 0
    good_coprime_pairs = 0

    # To be more efficient, we only need to iterate through q >= p
    # and multiply the result by 2 (and handle the p=q case).
    # This avoids counting (p,q) and (q,p) separately.
    
    for p in range(2, n + 1):
        for q in range(p, n + 1):
            if math.gcd(p, q) == 1:
                # This is a candidate pair.
                # If p != q, this represents two pairs: (p,q) and (q,p)
                count = 2 if p != q else 1
                total_coprime_pairs += count
                
                # Check if it's a "good" pair
                if p + q <= n + 1:
                    good_coprime_pairs += count
                    
    if total_coprime_pairs == 0:
        probability = 0
    else:
        probability = good_coprime_pairs / total_coprime_pairs

    print(f"The condition for a pair (p,q) to be good is: p + q <= {n+1}")
    print("\n--- Probability Calculation ---")
    print(f"Number of 'good' coprime pairs (p+q <= {n+1}): {good_coprime_pairs}")
    print(f"Total number of coprime pairs (2 <= p,q <= {n}): {total_coprime_pairs}")
    print(f"Pr({n}) = {good_coprime_pairs} / {total_coprime_pairs} = {probability:.6f}")
    print("\nAs n approaches infinity, this probability converges to 1/2.")

if __name__ == '__main__':
    # Using a reasonably large value for n for better approximation
    n_value = 2000
    calculate_probability(n_value)