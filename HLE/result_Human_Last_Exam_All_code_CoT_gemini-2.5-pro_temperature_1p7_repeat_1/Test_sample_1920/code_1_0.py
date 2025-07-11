import math
from sympy import factorint

def get_prime_factors(n):
    """Return a list of prime factors of n."""
    return list(factorint(n).keys())

def count_primitive_chars_order_dividing_k(k, primes):
    """
    Calculates the number of primitive Dirichlet characters of a square-free
    conductor d=p1*p2*... whose order divides k.
    """
    count_per_prime = 1
    # For a prime p, the number of primitive characters mod p with order dividing k is
    # gcd(k, p-1) - 1.
    # For all prime factors of d=53599, p-1 is divisible by k for any k|6.
    # So gcd(k, p-1) simplifies to k.
    if k == 0:
        return 0
    count_per_prime = k - 1
    
    # The total number is the product of counts for each prime factor.
    return count_per_prime ** len(primes)

def solve():
    """
    Finds the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6
    
    try:
        primes = get_prime_factors(d)
    except NameError:
        print("Please install sympy ('pip install sympy') for prime factorization.")
        # Fallback for the specific problem if sympy is not available
        primes = [7, 13, 19, 31]
    
    # Using inclusion-exclusion principle for order g=6.
    # Divisors of 6 are 1, 2, 3, 6.
    # Mobius values: mu(6/6)=mu(1)=1, mu(6/3)=mu(2)=-1, mu(6/2)=mu(3)=-1, mu(6/1)=mu(6)=1
    # Result = 1*N_6 - 1*N_3 - 1*N_2 + 1*N_1
    
    n_6 = count_primitive_chars_order_dividing_k(6, primes)
    n_3 = count_primitive_chars_order_dividing_k(3, primes)
    n_2 = count_primitive_chars_order_dividing_k(2, primes)
    n_1 = count_primitive_chars_order_dividing_k(1, primes)

    result = n_6 - n_3 - n_2 + n_1

    print(f"Let d = {d} and the order be g = {g}.")
    print(f"The prime factorization of d is: {primes}")
    print("\nThe number of primitive characters with order exactly 6 can be found using the inclusion-exclusion principle.")
    print("Let N_k be the number of primitive characters with conductor d whose order divides k.")
    print(f"N_6 = (6-1)^{len(primes)} = {n_6}")
    print(f"N_3 = (3-1)^{len(primes)} = {n_3}")
    print(f"N_2 = (2-1)^{len(primes)} = {n_2}")
    print(f"N_1 = (1-1)^{len(primes)} = {n_1}")
    
    print("\nThe total number is N_6 - N_3 - N_2 + N_1")
    print(f"= {n_6} - {n_3} - {n_2} + {n_1}")
    print(f"= {result}")

solve()