import math
from collections import defaultdict
import itertools

def get_prime_factorization_exponents(n, primes):
    """Calculates the exponents of specified primes in the factorization of n."""
    exponents = {p: 0 for p in primes}
    if n == 0: return exponents
    for p in primes:
        while n > 0 and n % p == 0:
            exponents[p] += 1
            n //= p
    return exponents

def check_n27():
    """
    Checks if it's possible to form a valid set for N=27.
    A key property is that for a prime p, if floor(N/p) < 4, no number in the grid
    can be a multiple of p. For N=27, this means p must be <= 6.75, so any
    number in the grid must only have prime factors {2, 3, 5}.
    """
    print("Step 1: Checking if N=27 is possible.")
    primes = {2, 3, 5}
    pool = []
    for i in range(1, 28):
        is_5_smooth = True
        factors = get_prime_factorization_exponents(i, {p for p in range(2, i + 1) if p > 5})
        if any(v > 0 for v in factors.values()):
            is_5_smooth = False
        if is_5_smooth:
            pool.append(i)
    
    print(f"The only allowed numbers for N=27 must be 5-smooth. Pool of numbers: {pool}")
    print(f"Size of pool: {len(pool)}. We need to select 16 numbers.")
    
    # We must select 16 numbers from this pool of 17.
    # Let's find the total exponent vector for the product of all 17 numbers.
    total_exponents = {p: 0 for p in primes}
    for num in pool:
        exponents = get_prime_factorization_exponents(num, primes)
        for p, e in exponents.items():
            total_exponents[p] += e
            
    total_exponents_mod4 = {p: e % 4 for p, e in total_exponents.items()}
    print(f"To form a valid set of 16, we must remove one number 'r' such that the")
    print(f"exponent vector of 'r' (mod 4) matches {total_exponents_mod4}.")

    found_candidate = False
    for num_to_remove in pool:
        exp_to_remove = get_prime_factorization_exponents(num_to_remove, primes)
        exp_to_remove_mod4 = {p: e % 4 for p, e in exp_to_remove.items()}
        if exp_to_remove_mod4 == total_exponents_mod4:
            found_candidate = True
            break
            
    if not found_candidate:
        print("No such number found to remove. Therefore, no subset of 16 has a product that is a perfect fourth power.")
        print("Conclusion: N=27 is not possible.\n")
        return False
    return True

def construct_n28():
    """
    Constructs a valid set and table for N=28.
    """
    print("Step 2: Proving N=28 is possible by construction.")
    s1 = [1, 2, 3, 4]
    s2 = [1, 5, 6, 7]
    
    numbers = sorted([i * j for i in s1 for j in s2])
    N = max(numbers)

    print(f"Using construction sets S1={s1} and S2={s2}.")
    print(f"The 16 distinct numbers are: {numbers}")
    print(f"The maximum value is N = {N}, so N=28 is achievable.")

    # A known construction for a multiplicative magic square of order 4
    # from sets s1 and s2.
    a1, a2, a3, a4 = s1
    b1, b2, b3, b4 = s2
    
    # Using a pair of orthogonal Latin squares
    table = [
        [a1*b1, a2*b2, a3*b3, a4*b4],
        [a2*b3, a1*b4, a4*b1, a3*b2],
        [a3*b4, a4*b3, a1*b2, a2*b1],
        [a4*b2, a3*b1, a2*b4, a1*b3]
    ]

    P = math.prod(s1) * math.prod(s2)
    print("\nA valid 4x4 table is:")
    for row in table:
        print(" ".join(f"{num:2d}" for num in row))

    print(f"\nThe product P of each row and column is: {P}")
    # Verification
    valid = True
    for i in range(4):
        row_prod = math.prod(table[i])
        if row_prod != P: valid = False
        col_prod = math.prod(table[j][i] for j in range(4))
        if col_prod != P: valid = False
    print(f"Verification successful: {valid}")


if __name__ == "__main__":
    if not check_n27():
        construct_n28()
        print("\nFinal Answer: The smallest N is 28.")
        print("<<<28>>>")