import math

def prime_factorize(n):
    """
    Returns a list of prime factors of a natural number.
    """
    factors = []
    # Check for divisibility by 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2
    if n > 1:
        factors.append(n)
    return factors

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    """
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

def solve():
    """
    Finds the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    order = 6
    print(f"Goal: Find the number of primitive Dirichlet characters of conductor d = {d} and order {order}.")
    print("-" * 20)

    # Step 1: Factorize the conductor
    print("Step 1: Factorize the conductor.")
    prime_factors = prime_factorize(d)
    factors_str = ' * '.join(map(str, prime_factors))
    print(f"d = {d} = {factors_str}")
    
    # Check if d is square-free, which simplifies the problem
    if len(prime_factors) != len(set(prime_factors)):
        print("Conductor is not square-free. This calculation is more complex.")
        return

    print("The conductor is square-free. A character is primitive if its components are non-trivial.")
    r = len(prime_factors)
    print(f"The number of prime factors is r = {r}.")
    
    # Step 2: Check conditions on prime factors
    print("\nStep 2: Verify conditions for characters of order 6.")
    print("For a character of order 6 to exist for a prime component p, 6 must divide p-1.")
    for p in prime_factors:
        print(f"For p = {p}, p-1 = {p-1}. (p-1) % 6 = {(p-1)%6}.")
    print("The condition holds for all prime factors.")

    # Step 3: Count characters using inclusion-exclusion
    print("\nStep 3: Count characters using the Principle of Inclusion-Exclusion.")
    
    # Total primitive characters whose order divides 6
    # For each prime p, the number of non-trivial characters with order dividing 6 is phi(2)+phi(3)+phi(6)
    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)
    choices_div_6 = phi_2 + phi_3 + phi_6
    total_div_6 = choices_div_6 ** r
    print(f"Total primitive characters with order dividing 6: ({phi_2} + {phi_3} + {phi_6})^{r} = {choices_div_6}^{r} = {total_div_6}")
    
    # Total primitive characters whose order divides 3
    # This means order of each component is 3 (since non-trivial)
    choices_div_3 = phi_3
    total_div_3 = choices_div_3 ** r
    print(f"Primitive characters with order dividing 3: {choices_div_3}^{r} = {total_div_3}")

    # Total primitive characters whose order divides 2
    # This means order of each component is 2
    choices_div_2 = phi_2
    total_div_2 = choices_div_2 ** r
    print(f"Primitive characters with order dividing 2: {choices_div_2}^{r} = {total_div_2}")
    
    # Characters with order dividing gcd(2,3)=1 is 0 as components are non-trivial.

    # Final calculation
    final_count = total_div_6 - total_div_3 - total_div_2
    
    print("\nFinal Calculation:")
    print("Number of characters of order 6 = (order divides 6) - (order divides 3) - (order divides 2)")
    print(f"{total_div_6} - {total_div_3} - {total_div_2} = {final_count}")

solve()