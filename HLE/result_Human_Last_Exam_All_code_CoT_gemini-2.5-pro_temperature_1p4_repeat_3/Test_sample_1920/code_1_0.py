import math

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
    d = 53599
    order_k = 6
    
    print(f"We want to find the number of primitive Dirichlet characters of conductor d = {d} and order {order_k}.")
    
    # Step 1: Factor the conductor d
    prime_factors_d = get_prime_factorization(d)
    primes = list(prime_factors_d.keys())
    num_prime_factors = len(primes)
    
    print(f"\nStep 1: Find the prime factorization of d.")
    print(f"The prime factorization of {d} is {' * '.join(map(str, primes))}.")
    
    is_square_free = all(v == 1 for v in prime_factors_d.values())
    if not is_square_free:
        print("The conductor d is not square-free. The calculation is more complex and not covered here.")
        return
        
    print("Since d is square-free, a character modulo d is primitive if and only if its component for each prime factor is non-trivial (i.e., has order > 1).")

    print(f"\nStep 2: Understand the structure of the character group.")
    print(f"The order of a character chi = (chi_p1, ..., chi_p{num_prime_factors}) is lcm(ord(chi_p1), ..., ord(chi_p{num_prime_factors})).")
    print(f"We need this lcm to be {order_k}.")
    print("This implies that the order of each component must be a divisor of 6, i.e., in {1, 2, 3, 6}.")
    print("The primitivity condition means the order of each component must be > 1, so the possible orders for each component are {2, 3, 6}.")
    
    print("\nStep 3: Use the Principle of Inclusion-Exclusion.")
    print("Let N(m) be the number of primitive characters whose order divides m.")
    print(f"The number of characters of order exactly {order_k} is: N(6) - N(3) - N(2) + N(1).")

    # Step 4: Calculate each term N(m)
    print("\nStep 4: Calculate N(m) for m = 6, 3, 2, 1.")
    
    num_s = {}
    
    # This calculation relies on the fact that for each prime p in the factorization of d, (p-1) is a multiple of 6.
    # 7-1=6, 13-1=12, 19-1=18, 31-1=30.
    for m in [6, 3, 2, 1]:
        # For a prime p, the number of characters mod p whose order divides m is m, provided m | (p-1).
        # Number of non-trivial characters is m-1.
        num_non_trivial_per_prime = m - 1
        term_result = num_non_trivial_per_prime ** num_prime_factors
        num_s[m] = term_result
        
        print(f"\nCalculating N({m}):")
        print(f"The number of characters mod p with order dividing {m} is {m}.")
        print(f"The number of non-trivial characters (order > 1) is {m} - 1 = {num_non_trivial_per_prime}.")
        print(f"Since there are {num_prime_factors} prime factors, N({m}) = {num_non_trivial_per_prime}^{num_prime_factors} = {term_result}.")

    # Final calculation
    result = num_s[6] - num_s[3] - num_s[2] + num_s[1]
    
    print("\nStep 5: Final Calculation.")
    print(f"The number of primitive Dirichlet characters of conductor {d} and order {order_k} is:")
    print(f"{num_s[6]} - {num_s[3]} - {num_s[2]} + {num_s[1]} = {result}")

if __name__ == "__main__":
    main()