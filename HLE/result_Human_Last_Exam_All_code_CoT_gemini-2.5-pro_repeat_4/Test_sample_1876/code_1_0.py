import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def check_set_product(number_set):
    """
    Checks if the product of numbers in a set is a perfect fourth power.
    """
    total_exponents = {}
    print(f"The set of 16 distinct integers is:\n{sorted(list(number_set))}\n")
    
    for num in number_set:
        factors = get_prime_factorization(num)
        for prime, exponent in factors.items():
            total_exponents[prime] = total_exponents.get(prime, 0) + exponent
            
    print("Sum of exponents for each prime in the product of the set:")
    is_perfect_fourth = True
    for prime, exponent in sorted(total_exponents.items()):
        print(f"Prime {prime}: exponent = {exponent}", end="")
        if exponent % 4 == 0:
            print(" (divisible by 4)")
        else:
            print(f" (NOT divisible by 4, remainder is {exponent % 4})")
            is_perfect_fourth = False
            
    print("\nConclusion:")
    if is_perfect_fourth:
        print("The product of all numbers in the set is a perfect fourth power.")
    else:
        print("The product of all numbers in the set is NOT a perfect fourth power.")

# The set of 16 numbers with N=35
final_set = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 16, 25, 27, 28, 35}

check_set_product(final_set)

# The question asks for N, the maximum value in the set.
N = max(final_set)
print(f"\nThe smallest N is {N}.")
