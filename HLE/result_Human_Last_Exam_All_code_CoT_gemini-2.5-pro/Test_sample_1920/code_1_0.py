import math

def get_prime_factorization(n):
    """
    Returns a list of prime factors of a natural number n.
    """
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.append(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    return factors

def phi(n):
    """
    Calculates Euler's totient function phi(n).
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
    return int(result)

def solve():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d={d} and order g={g}.\n")

    # Step 1: Prime factorization of d
    primes = get_prime_factorization(d)
    print(f"Step 1: The prime factorization of d={d} is {' * '.join(map(str, primes))}.")
    print("The conductor is square-free.\n")

    # Step 2: Define conditions for component characters
    print("Step 2: A primitive character chi modulo d is a product of non-principal characters chi_i modulo each prime factor p_i.")
    print(f"The order of chi, lcm(ord(chi_1), ord(chi_2), ord(chi_3), ord(chi_4)), must be {g}.")
    print(f"This means ord(chi_i) must be a divisor of {g} and greater than 1.")
    
    divisors_g = [k for k in range(1, g + 1) if g % k == 0]
    possible_orders = [k for k in divisors_g if k > 1]
    print(f"The possible orders for each component chi_i are {possible_orders}.\n")

    # Step 3: Count characters for each possible order
    print("Step 3: Count the number of characters for each possible order modulo the prime factors.")
    phi_values = {k: phi(k) for k in possible_orders}
    print("The number of characters of order k modulo a prime p is phi(k), provided k divides (p-1).")
    
    all_valid = True
    for p in primes:
        if (p - 1) % g != 0:
            print(f"Warning: For prime {p}, not all orders dividing {g} are possible since {g} does not divide {p-1}.")
            all_valid = False
            break
    if all_valid:
        print(f"For all prime factors p of {d}, (p-1) is divisible by {g}, so all orders {possible_orders} are possible.")
    
    for k, v in phi_values.items():
        print(f"  - Number of characters of order {k} is phi({k}) = {v}")
    print("")

    # Step 4: Use Inclusion-Exclusion Principle
    print("Step 4: Use the Principle of Inclusion-Exclusion to find the final count.")
    
    # Total primitive characters with order dividing 6
    num_choices_div6 = sum(phi_values.values())
    total_div6 = num_choices_div6 ** len(primes)
    
    # Total primitive characters with order dividing 3
    num_choices_div3 = phi_values.get(3, 0)
    total_div3 = num_choices_div3 ** len(primes)
    
    # Total primitive characters with order dividing 2
    num_choices_div2 = phi_values.get(2, 0)
    total_div2 = num_choices_div2 ** len(primes)
    
    # Total primitive characters with order dividing 1 (is 0, as order must be > 1)
    total_div1 = 0

    result = total_div6 - total_div3 - total_div2 + total_div1

    print("The number of characters of order exactly 6 is given by the formula:")
    print("N(6) = (Num. with order dividing 6) - (Num. with order dividing 3) - (Num. with order dividing 2) + (Num. with order dividing 1)")
    print(f"The number of choices for one component's order to be in {{2, 3, 6}} is {phi(2)} + {phi(3)} + {phi(6)} = {num_choices_div6}.")
    print(f"The number of choices for one component's order to be in {{3}} is {phi(3)} = {num_choices_div3}.")
    print(f"The number of choices for one component's order to be in {{2}} is {phi(2)} = {num_choices_div2}.")
    print(f"Since there are {len(primes)} prime factors, we have:")
    print(f"N(6) = {num_choices_div6}^{len(primes)} - {num_choices_div3}^{len(primes)} - {num_choices_div2}^{len(primes)} + 0")
    print(f"N(6) = {total_div6} - {total_div3} - {total_div2} = {result}")

    print(f"\nFinal Answer: The number of primitive Dirichlet characters is {result}.")

solve()