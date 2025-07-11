import math

def prime_factors(n):
    """
    Returns a list of prime factors of a positive integer.
    """
    factors = []
    temp_n = n
    d = 2
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.append(d)
            temp_n //= d
        else:
            d += 1
    if temp_n > 1:
       factors.append(temp_n)
    return factors

def phi(n):
    """
    Calculates Euler's totient function.
    """
    if n == 1:
        return 1
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

def solve():
    """
    Finds the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    k = 6

    print(f"The conductor is d = {d}.")
    print(f"The required order is k = {k}.")

    # Step 1: Prime factorization of the conductor d
    factors_d = prime_factors(d)
    r = len(factors_d)
    print(f"\nStep 1: The prime factorization of d is {d} = {' * '.join(map(str, factors_d))}.")
    print("A primitive character of conductor d is a product of r non-principal characters, one for each prime factor.")
    
    print("\nStep 2: We need to find the number of tuples of non-principal characters (chi_1, ..., chi_r) such that")
    print(f"lcm(ord(chi_1), ..., ord(chi_r)) = {k}.")
    print("This will be solved using the Principle of Inclusion-Exclusion.")
    print(f"The number of characters is N(6) - N(3) - N(2) + N(1), where N(m) is the number of combinations")
    print(f"of characters where the lcm of their orders divides m.")

    # Step 3: Calculate the terms for the formula
    print("\nStep 3: Calculating N(m) for m = 6, 3, 2, 1.")
    
    # N(6)
    num_choices_6 = sum(phi(j) for j in [2, 3, 6])
    term_N6 = num_choices_6 ** r
    print(f"The number of non-principal characters mod p_i whose order divides 6 is phi(2)+phi(3)+phi(6) = {phi(2)}+{phi(3)}+{phi(6)} = {num_choices_6}.")
    print(f"With r={r} prime factors, N(6) = {num_choices_6}^{r} = {term_N6}.")

    # N(3)
    num_choices_3 = phi(3)
    term_N3 = num_choices_3 ** r
    print(f"\nThe number of non-principal characters mod p_i whose order divides 3 is phi(3) = {phi(3)} = {num_choices_3}.")
    print(f"N(3) = {num_choices_3}^{r} = {term_N3}.")
    
    # N(2)
    num_choices_2 = phi(2)
    term_N2 = num_choices_2 ** r
    print(f"\nThe number of non-principal characters mod p_i whose order divides 2 is phi(2) = {phi(2)} = {num_choices_2}.")
    print(f"N(2) = {num_choices_2}^{r} = {term_N2}.")

    # N(1)
    term_N1 = 0
    print("\nThe number of non-principal characters whose order divides 1 is 0. So N(1) = 0.")

    # Step 4: Final Calculation
    result = term_N6 - term_N3 - term_N2 + term_N1
    print("\nStep 4: Final calculation using inclusion-exclusion.")
    print("Number of characters = N(6) - N(3) - N(2) + N(1)")
    print(f"= {term_N6} - {term_N3} - {term_N2} + {term_N1}")
    print(f"= {result}")

solve()