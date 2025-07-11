import math

def get_prime_factorization(n):
    """
    Finds the unique prime factors of an integer n.
    """
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
    return list(factors)

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    """
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return int(result)

def num_primitive_chars_prime_conductor(p, g):
    """
    Calculates the number of primitive characters of conductor p and order g.
    For a prime conductor p, this is phi(g) if g > 1 and g divides p-1.
    """
    if g == 1:
        return 0
    if (p - 1) % g == 0:
        return phi(g)
    return 0

def calculate_A(k, primes):
    """
    Calculates A(k), the number of primitive characters with conductor d
    whose order divides k.
    """
    total_prod = 1
    for p in primes:
        # Sum of primitive characters for a single prime factor p,
        # for orders g that divide k.
        current_sum = 0
        # Iterate through all divisors of k
        for g in range(1, k + 1):
            if k % g == 0:
                current_sum += num_primitive_chars_prime_conductor(p, g)
        total_prod *= current_sum
    return total_prod

def main():
    """
    Main function to solve the problem.
    """
    d = 53599
    order = 6

    # Step 1: Factorize the conductor
    primes = get_prime_factorization(d)

    # Step 4 & 5: Calculate A(k) for divisors of the order
    # The divisors of 6 are 1, 2, 3, 6.
    A6 = calculate_A(6, primes)
    A3 = calculate_A(3, primes)
    A2 = calculate_A(2, primes)
    A1 = calculate_A(1, primes)

    # Step 6: Final calculation using inclusion-exclusion
    result = A6 - A3 - A2 + A1

    print(f"The conductor is d = {d}.")
    print(f"The prime factorization of d is: {' * '.join(map(str, sorted(primes)))}")
    print("\nTo find the number of characters of order exactly 6, we use the inclusion-exclusion principle.")
    print("Let A(k) be the number of characters whose order divides k.")
    print("The total number of characters of order 6 is A(6) - A(3) - A(2) + A(1).\n")
    print(f"Number of characters with order dividing 6, A(6) = {A6}")
    print(f"Number of characters with order dividing 3, A(3) = {A3}")
    print(f"Number of characters with order dividing 2, A(2) = {A2}")
    print(f"Number of characters with order dividing 1, A(1) = {A1}\n")
    print("The final number of primitive Dirichlet characters is calculated as:")
    print(f"{A6} - {A3} - {A2} + {A1} = {result}")

if __name__ == "__main__":
    main()
