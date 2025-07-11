import math

def count_allowed_pairs():
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a, b <= 1000.

    The conditions for a pair (a,b) to be allowed are:
    1. a = 1 or b = 1.
    2. a is a prime p, and b = p^k for k >= 1.
    3. b is a prime p, and a = p^k for k >= 1.

    This function counts the pairs satisfying these conditions using the
    Principle of Inclusion-Exclusion.
    """
    N = 1000

    # Let C1, C2, C3 be the sets of pairs for each condition.
    # Total = |C1| + |C2| + |C3| - (|C1 n C2| + |C1 n C3| + |C2 n C3|) + |C1 n C2 n C3|

    # Count for Condition 1: a = 1 or b = 1
    # Pairs with a=1: 1000. Pairs with b=1: 1000. The pair (1,1) is in both sets.
    count1 = 1000 + 1000 - 1

    # For conditions 2 and 3, we need primes up to N.
    # Generate primes using a Sieve of Eratosthenes.
    is_prime = [True] * (N + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(N)) + 1):
        if is_prime[i]:
            for multiple in range(i * i, N + 1, i):
                is_prime[multiple] = False
    primes = [i for i, p_bool in enumerate(is_prime) if p_bool]
    num_primes = len(primes)

    # Count for Condition 2: a = p, b = p^k
    # where p is a prime and k >= 1, with a <= N and b <= N.
    count2 = 0
    for p in primes:
        # a = p is always <= N
        power_val = p
        while power_val <= N:
            count2 += 1
            # Check for potential overflow before multiplication
            if N // p < power_val:
                break
            power_val *= p

    # Count for Condition 3: b = p, a = p^k
    # This is symmetric to condition 2.
    count3 = count2
    
    # --- Intersection Counts ---

    # Intersection of C1 and C2 (or C1 and C3) is empty.
    # For C2 and C3, a > 1 and b > 1.
    count12 = 0
    count13 = 0

    # Intersection of C2 and C3: (a,b) is of the form (p, p^k) AND (q^m, q).
    # This implies a=p=q^m and b=p^k=q.
    # a=p, p is prime. a=q^m, q is prime. This forces m=1 and p=q.
    # b=q, q is prime. b=p^k, p is prime. This forces k=1 and p=q.
    # The intersection consists of pairs (p,p) where p is a prime <= N.
    count23 = num_primes

    # Intersection of C1, C2, C3 is empty.
    count123 = 0

    # Applying the Principle of Inclusion-Exclusion
    total_allowed = count1 + count2 + count3 - (count12 + count13 + count23) + count123
    
    print(f"Number of pairs from Condition 1 (|C1|): {count1}")
    print(f"Number of pairs from Condition 2 (|C2|): {count2}")
    print(f"Number of pairs from Condition 3 (|C3|): {count3}")
    print(f"Number of pairs in C2 and C3 (|C2 n C3|): {count23}")
    print("\nFinal calculation:")
    print(f"Total = |C1| + |C2| + |C3| - |C2 n C3|")
    print(f"Total = {count1} + {count2} + {count3} - {count23}")
    print(f"Total = {count1 + count2 + count3 - count23}")

    print(f"\nThere are {total_allowed} allowed pairs.")
    print(f"<<<{total_allowed}>>>")

count_allowed_pairs()