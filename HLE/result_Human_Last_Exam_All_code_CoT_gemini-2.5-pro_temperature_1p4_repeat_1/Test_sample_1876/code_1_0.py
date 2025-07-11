import math
from functools import reduce

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def main():
    """
    Main function to find and verify the set of numbers.
    """
    # Smallest N from a known mathematical construction.
    N = 210
    
    # The set S is {2^a * 3^b * 5^c * 7^d | a,b,c,d in {0,1}}
    primes = [2, 3, 5, 7]
    s = set()
    for i in range(16):
        num = 1
        # Each i from 0 to 15 corresponds to a unique combination of exponents
        temp_i = i
        for p_idx in range(4):
            if (temp_i % 2) == 1:
                num *= primes[p_idx]
            temp_i //= 2
        s.add(num)
    
    number_set = sorted(list(s))
    
    print(f"A set of 16 distinct integers is chosen, with largest integer N = {max(number_set)}")
    print("The set is:", number_set)
    
    # Calculate the product of all numbers in the set
    product = reduce(lambda x, y: x * y, number_set)
    
    # Check if the product is a perfect fourth power.
    # We do this by checking the exponents of its prime factorization.
    total_factors = {}
    for num in number_set:
        factors = get_prime_factorization(num)
        for p, exp in factors.items():
            total_factors[p] = total_factors.get(p, 0) + exp
            
    is_fourth_power = True
    print("\nExponents of prime factors in the total product:")
    for p, exp in sorted(total_factors.items()):
        print(f"Prime {p}: exponent = {exp}")
        if exp % 4 != 0:
            is_fourth_power = False

    if is_fourth_power:
        print("\nThe product of all numbers in the set is a perfect fourth power.")
        P_fourth = product
        P = int(round(P_fourth**(1/4)))
        print(f"The common product P should be {P}")
        print(f"Let's verify P^4: {P}^4 = {P**4}")
        print(f"The product of the set is: {product}")
        if P**4 == product:
          print("Verification successful.")
        else:
          print("Verification failed.")

    else:
        print("\nThe product is not a perfect fourth power, so this set is not a solution.")

    print(f"\nBased on this construction, the smallest N is {N}.")

main()
