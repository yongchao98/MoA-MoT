import math

def get_prime_factorization(n):
    """Computes the prime factorization of n."""
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

def get_divisors(n):
    """Gets all divisors of n from its prime factorization."""
    factors = list(get_prime_factorization(n).items())
    divs = [1]
    for p, a in factors:
        # For each prime factor p with exponent a, extend the list of divisors
        divs.extend([d * (p**i) for i in range(1, a + 1) for d in divs])
    return sorted(list(set(divs)))

def phi(n):
    """Calculates Euler's totient function phi(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return result

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b) if a != 0 and b != 0 else 0

def solve_character_count():
    """
    Finds the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.\n")

    # Step 1: Factorize the conductor d
    d_factors = get_prime_factorization(d)
    if len(d_factors) != 2 or any(v != 1 for v in d_factors.values()):
        print(f"Error: This method is for conductors that are a product of two distinct primes.")
        return

    p1, p2 = d_factors.keys()
    print(f"Step 1: The conductor d = {d} is a product of two distinct primes: p1 = {p1}, p2 = {p2}.")

    print("\nA primitive character chi of conductor d = p1*p2 is a product chi = chi_1 * chi_2,")
    print("where chi_1 and chi_2 are non-trivial characters modulo p1 and p2 respectively.")
    print("This implies their orders g1 and g2 must both be greater than 1.\n")

    # Step 2: Set up conditions on the orders g1 and g2
    phi1 = p1 - 1
    phi2 = p2 - 1
    print(f"Step 2: The orders must satisfy the following conditions:")
    print(f"  - g1 > 1 and g2 > 1")
    print(f"  - g1 must divide (p1 - 1) = {phi1}")
    print(f"  - g2 must divide (p2 - 1) = {phi2}")
    print(f"  - lcm(g1, g2) = {g}\n")

    # Step 3: Find all pairs of orders (g1, g2) satisfying the conditions
    divs1 = get_divisors(phi1)
    divs2 = get_divisors(phi2)
    
    valid_pairs = []
    for g1 in divs1:
        if g1 <= 1:
            continue
        for g2 in divs2:
            if g2 <= 1:
                continue
            if lcm(g1, g2) == g:
                valid_pairs.append((g1, g2))

    print(f"Step 3: The valid pairs of orders (g1, g2) found are: {valid_pairs}\n")

    # Step 4: Calculate the total number of characters
    total_count = 0
    full_equation_parts = []
    if not valid_pairs:
        print("No valid pairs found, so the total number of characters is 0.")
    else:
        print("Step 4: For each pair (g1, g2), the number of characters is phi(g1) * phi(g2).")
        for g1, g2 in valid_pairs:
            num_chi1 = phi(g1)
            num_chi2 = phi(g2)
            count = num_chi1 * num_chi2
            total_count += count
            
            print(f"\n- For (g1, g2) = ({g1}, {g2}):")
            print(f"  Number of characters of order {g1} mod {p1} = phi({g1}) = {num_chi1}")
            print(f"  Number of characters of order {g2} mod {p2} = phi({g2}) = {num_chi2}")
            print(f"  Contribution to total = {num_chi1} * {num_chi2} = {count}")
            full_equation_parts.append(f"{num_chi1} * {num_chi2}")
        
        final_equation = " + ".join(full_equation_parts)
        print("\nFinal Calculation:")
        print(f"Total number = {final_equation} = {total_count}")

solve_character_count()