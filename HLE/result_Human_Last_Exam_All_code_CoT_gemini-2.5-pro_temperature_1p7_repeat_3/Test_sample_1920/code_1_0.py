import math

def prime_factorize(n):
    """Returns a list of prime factors of n."""
    factors = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors.append(d)
            temp_n //= d
        d += 1
    if temp_n > 1:
       factors.append(temp_n)
    return factors

def phi(n):
    """Computes Euler's totient function."""
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
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {k}.\n")

    # Step 1: Factorize the conductor
    prime_factors = prime_factorize(d)
    r = len(prime_factors)
    print(f"Step 1: The conductor d = {d} is a product of r = {r} distinct primes:")
    print(f"d = {' * '.join(map(str, prime_factors))}\n")
    print("For the character to be primitive, it must be a product of non-trivial characters modulo each prime factor.\n")
    
    # Step 2: Set up the counting based on order k=6
    print("Step 2: The order of the character is the lcm of the orders of its components.")
    print("For the order to be 6, the component orders must divide 6, i.e., belong to {2, 3, 6} (since they must be > 1).\n")

    # Step 3: Use inclusion-exclusion
    # First, count characters where component orders are in {2, 3, 6}
    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)
    
    print("Step 3: Count the number of choices for each component character.")
    print(f"Number of characters of order 2: phi(2) = {phi_2}")
    print(f"Number of characters of order 3: phi(3) = {phi_3}")
    print(f"Number of characters of order 6: phi(6) = {phi_6}\n")

    # Total combinations where order divides 6 (and > 1)
    choices_per_comp_div_6 = phi_2 + phi_3 + phi_6
    total_div_6 = choices_per_comp_div_6 ** r
    
    # Combinations where lcm divides 3 (all component orders must be 3)
    choices_per_comp_div_3 = phi_3
    total_div_3 = choices_per_comp_div_3 ** r

    # Combinations where lcm divides 2 (all component orders must be 2)
    choices_per_comp_div_2 = phi_2
    total_div_2 = choices_per_comp_div_2 ** r

    print("Step 4: Apply the Principle of Inclusion-Exclusion.")
    print("Total combinations with component orders dividing 6:")
    print(f"({phi_2} + {phi_3} + {phi_6})^{r} = {choices_per_comp_div_6}^{r} = {total_div_6}\n")

    print("Subtract combinations where the lcm is a proper divisor of 6.")
    print("Combinations with lcm dividing 3 (orders must be 3):")
    print(f"({phi_3})^{r} = {choices_per_comp_div_3}^{r} = {total_div_3}")
    print("Combinations with lcm dividing 2 (orders must be 2):")
    print(f"({phi_2})^{r} = {choices_per_comp_div_2}^{r} = {total_div_2}\n")

    # Final result
    result = total_div_6 - total_div_3 - total_div_2

    print("Final Result:")
    print("The final count is obtained by the formula:")
    print(f"N = ({phi_2} + {phi_3} + {phi_6})^{r} - ({phi_3})^{r} - ({phi_2})^{r}")
    print(f"N = ({choices_per_comp_div_6})^{r} - ({choices_per_comp_div_3})^{r} - ({choices_per_comp_div_2})^{r}")
    print(f"N = {total_div_6} - {total_div_3} - {total_div_2}")
    print(f"N = {result}")

solve()