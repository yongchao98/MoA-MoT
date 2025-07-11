import math

def get_prime_factorization(n):
    """
    Returns a list of unique prime factors of n.
    Raises ValueError if n is not square-free.
    """
    factors = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.append(d)
            temp_n //= d
            if temp_n % d == 0:
                raise ValueError(f"{n} is not square-free, this method is for square-free conductors.")
        d += 1
    if temp_n > 1:
        factors.append(temp_n)
    return factors

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

def solve():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {k}.\n")

    try:
        # Step 1: Factorize the conductor d
        prime_factors = get_prime_factorization(d)
        r = len(prime_factors)
        print(f"Step 1: The conductor d = {d} is square-free.")
        print(f"The prime factorization is d = {' * '.join(map(str, prime_factors))}.")
        print(f"The number of prime factors is r = {r}.\n")

        # Step 2: Verify conditions for existence of characters of required order
        print(f"Step 2: For a character component of order m to exist modulo a prime p, m must divide p-1.")
        print(f"We need characters of order dividing 6, so we check if 6 divides p-1 for each prime factor p.")
        for p in prime_factors:
            if (p - 1) % k != 0:
                print(f"Error: For prime factor {p}, {p-1} is not divisible by {k}.")
                print("Cannot form characters of order 6. The number of such characters is 0.")
                return
        print(f"Condition satisfied: For all prime factors p, p-1 is divisible by {k}.\n")

        # Step 3: Count characters using inclusion-exclusion
        print("Step 3: We count the number of combinations of characters using inclusion-exclusion.")
        
        # Count non-principal characters whose order divides 6
        divisors_6 = [2, 3, 6]
        count_div_6 = sum(phi(j) for j in divisors_6)
        total_div_6 = count_div_6 ** r
        print(f"Number of non-principal characters mod p with order dividing 6 is phi(2)+phi(3)+phi(6) = {phi(2)}+{phi(3)}+{phi(6)} = {count_div_6}.")
        print(f"Total primitive characters with order dividing 6 is {count_div_6}^{r} = {total_div_6}.\n")

        # Count non-principal characters whose order divides 2 (i.e., is 2)
        divisors_2 = [2]
        count_div_2 = sum(phi(j) for j in divisors_2)
        total_div_2 = count_div_2 ** r
        print(f"Number of characters mod p with order 2 is phi(2) = {phi(2)}.")
        print(f"Total primitive characters with order dividing 2 is {count_div_2}^{r} = {total_div_2}.\n")

        # Count non-principal characters whose order divides 3 (i.e., is 3)
        divisors_3 = [3]
        count_div_3 = sum(phi(j) for j in divisors_3)
        total_div_3 = count_div_3 ** r
        print(f"Number of characters mod p with order 3 is phi(3) = {phi(3)}.")
        print(f"Total primitive characters with order dividing 3 is {count_div_3}^{r} = {total_div_3}.\n")

        # Final Calculation
        result = total_div_6 - total_div_2 - total_div_3
        print("Step 4: The number of characters of order 6 is the total number of characters with order dividing 6,")
        print("minus those with order dividing 2, minus those with order dividing 3.")
        
        # The problem asks to output the numbers in the final equation.
        print("\nFinal Calculation:")
        print(f"{count_div_6}^{r} - {count_div_2}^{r} - {count_div_3}^{r} = {total_div_6} - {total_div_2} - {total_div_3} = {result}")

    except ValueError as e:
        print(e)

solve()