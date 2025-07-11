import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
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

def get_divisors(n):
    """Returns a list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve():
    """Main function to solve the problem."""
    d = 53599
    order_k = 6

    print(f"We want to find the number of primitive Dirichlet characters of conductor d = {d} and order {order_k}.")
    print("-" * 20)

    # Step 1: Prime factorization
    print("Step 1: Find the prime factorization of d.")
    factors_d = get_prime_factorization(d)
    primes = sorted(factors_d.keys())
    print(f"The prime factorization of {d} is {' * '.join(map(str, primes))}.")
    print("Since d is square-free, a character is primitive with conductor d if and only if it is a product of non-principal characters for each prime factor.\n")
    
    # Step 2: Explain method
    print("Step 2: Use Inclusion-Exclusion to find characters of order exactly 6.")
    print("The number of such characters is given by the formula: N(6) - N(3) - N(2) + N(1), where N(m) is the number of primitive character tuples whose order's lcm divides m.\n")
    
    # Step 3: Calculate N(m) for divisors of k
    print("Step 3: Calculate N(m) for each divisor m of 6.")
    
    divs_of_k = get_divisors(order_k)
    phi_vals = {i: phi(i) for i in range(1, order_k + 1)}
    
    term_values = {}

    for m in divs_of_k:
        print(f"--- Calculating N({m}) ---")
        
        counts_per_prime = []
        for p in primes:
            p_minus_1 = p - 1
            common_divisors_gcd = math.gcd(m, p_minus_1)
            
            count = 0
            # Iterate through divisors of gcd(m, p-1)
            divs_of_gcd = get_divisors(common_divisors_gcd)
            for k_div in divs_of_gcd:
                if k_div > 1:  # Must be primitive (order > 1)
                    count += phi_vals.get(k_div, 0)
            
            counts_per_prime.append(count)
            print(f"For p={p}, the number of primitive characters with order dividing {m} is {count}.")
        
        N_m = 1
        for count in counts_per_prime:
            N_m *= count
        
        term_values[m] = N_m
        print(f"N({m}) = {' * '.join(map(str, counts_per_prime))} = {N_m}\n")

    # Step 4: Final calculation
    result = term_values.get(6, 0) - term_values.get(3, 0) - term_values.get(2, 0) + term_values.get(1, 0)
    
    print("Step 4: Final Calculation.")
    print(f"The total number of primitive characters of order 6 is:")
    print(f"N(6) - N(3) - N(2) + N(1) = {term_values[6]} - {term_values[3]} - {term_values[2]} + {term_values[1]}")
    print(f"= {result}")

solve()