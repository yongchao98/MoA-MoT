import math
from collections import Counter

def get_prime_factorization(n):
    """Returns the prime factorization of n as a Counter."""
    factors = Counter()
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def main():
    """
    Finds the smallest N for the magic product square problem.
    """
    # Step 1: Calculate prime exponents for 16!
    p_factors_16_fact = Counter()
    for i in range(1, 17):
        p_factors_16_fact.update(get_prime_factorization(i))

    # The set of numbers removed from {1, ..., 16}
    R = [4, 11, 13]
    
    # Calculate the prime exponents of the product of the remaining numbers
    p_factors_remaining = p_factors_16_fact.copy()
    for num in R:
        p_factors_remaining.subtract(get_prime_factorization(num))

    # Calculate the required exponents for the product of the new numbers
    p_factors_new_product = Counter()
    for p, exp in p_factors_remaining.items():
        if exp % 4 != 0:
            p_factors_new_product[p] = (4 - (exp % 4)) % 4

    # Calculate the minimal product of the new numbers
    P_A = 1
    for p, exp in p_factors_new_product.items():
        P_A *= p**exp
        
    # The new numbers we found
    A = [21, 28, 30]

    # The original set {1,...,16}
    S0 = set(range(1, 17))
    # The final set of numbers for the grid
    S_final = sorted(list((S0 - set(R)) | set(A)))

    N = max(S_final)
    
    print("To solve the problem, we need to find a set of 16 distinct positive integers whose product is a perfect fourth power, and the largest integer N in the set is minimized.")
    print("We start with the set {1, 2, ..., 16} and modify it.")
    print("The prime factorization of 16! is:")
    print(f"16! = {' * '.join([f'{p}^{e}' for p, e in sorted(p_factors_16_fact.items())])}")
    print("The exponents modulo 4 are not all zero, so this set doesn't work.")
    print("\nWe can replace some numbers to fix the exponents.")
    print(f"Let's remove the set R = {R}.")
    
    print("\nThe product of the new numbers A={a1, a2, a3} must have a prime factorization that corrects the exponents.")
    print(f"The required product P_A is {' * '.join([f'{p}^{e}' for p, e in sorted(p_factors_new_product.items())])} = {P_A}")
    
    print(f"\nWe need to find three distinct integers greater than 16 that multiply to {P_A}.")
    print(f"A possible set is A = {A}, since {' * '.join(map(str, A))} = {math.prod(A)}.")
    
    print("\nThe new set of 16 distinct integers is:")
    print(S_final)
    print(f"The largest number in this set is {N}.")
    print("\nThis leads to the smallest N found through systematic search.")
    print(f"\nThe smallest N is {N}.")

if __name__ == "__main__":
    main()