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

def get_divisors(n):
    """
    Returns a sorted list of all divisors of n.
    """
    if n == 0:
        return []
    divs = {1}
    factors = get_prime_factorization(n)
    primes = list(factors.keys())
    
    def generate(k, current_div):
        if k == len(primes):
            divs.add(current_div)
            return
        
        p = primes[k]
        for i in range(factors[p] + 1):
            generate(k + 1, current_div * (p**i))

    generate(0, 1)
    return sorted(list(divs))

def mobius_mu(n):
    """
    Calculates the Mobius mu function for n.
    """
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def main():
    """
    Calculates the number of primitive Dirichlet characters for a given
    conductor d and order g.
    """
    d = 53599
    g = 6

    # Step 1: Find prime factorization of d
    d_factors = get_prime_factorization(d)
    primes = list(d_factors.keys())
    
    # Check if d is square-free
    is_square_free = all(exp == 1 for exp in d_factors.values())
    if not is_square_free:
        print(f"Conductor d={d} is not square-free. This method is for square-free conductors.")
        return

    print(f"The conductor is d = 53599.")
    print(f"Its prime factorization is: {' * '.join(map(str, primes))}")
    print(f"The desired order is g = {g}.")
    print("-" * 30)
    print("We calculate S(j) for each divisor j of 6.")
    print("S(j) is the number of character tuples with order dividing j.")
    print("-" * 30)

    # Function to calculate S(j)
    def calculate_S(j, prime_list):
        term = 1
        for p in prime_list:
            term *= (math.gcd(j, p - 1) - 1)
        return term

    g_divisors = get_divisors(g)
    total_count = 0
    
    # We will build the equation string
    equation_parts = []
    
    # Store S(j) values for printing
    s_values = {}

    for j in g_divisors:
        s_j = calculate_S(j, primes)
        s_values[j] = s_j
        mu_val = mobius_mu(g // j)
        
        if mu_val == 0:
            continue
            
        term = mu_val * s_j
        total_count += term
        
        equation_parts.append(f"({mu_val}) * {s_j}")

    # Print S(j) values
    for j in sorted(s_values.keys()):
        print(f"S({j}) = {s_values[j]}")
    
    print("-" * 30)
    print("The total number of primitive characters of order 6 is calculated using the formula:")
    print(f"N = sum_{{j|{g}}} mu({g}/j) * S(j)")
    
    final_equation = " + ".join(equation_parts)
    print(f"N = {final_equation} = {total_count}")

main()
<<<608>>>