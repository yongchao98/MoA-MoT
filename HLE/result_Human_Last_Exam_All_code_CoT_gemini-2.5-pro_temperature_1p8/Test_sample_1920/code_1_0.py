import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def get_divisors(n):
    """Returns a list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def phi(n):
    """Computes Euler's totient function phi(n)."""
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
    Calculates the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    d_factors = get_prime_factorization(d)
    primes_of_d = list(d_factors.keys())
    
    # A character is primitive for a square-free conductor d if it's a product of 
    # primitive characters for each prime factor of d.
    # The conductor d must be square-free.
    for p, a in d_factors.items():
        if a > 1:
            print(f"Conductor d={d} is not square-free. This method does not apply.")
            return

    # Using Principle of Inclusion-Exclusion for the order g=6.
    # The divisors of 6 are 1, 2, 3, 6.
    # N(6) = N(orders divide 6) - N(orders divide 3) - N(orders divide 2) + N(orders divide 1)
    
    divisors_of_g = get_divisors(g)
    # The inclusion-exclusion is over the prime factors of g, not all divisors.
    # lcm(orders)=6 means (for all p, order(p)|6) AND (not for all p, order(p)|3) AND (not for all p, order(p)|2)
    # which is N(lcm|6) - N(lcm|3) - N(lcm|2) + N(lcm|gcd(2,3)=1)
    
    terms = {}

    for k in divisors_of_g:
        total_prod = 1
        num_choices_per_prime = []
        for p in primes_of_d:
            # Count primitive characters for conductor p whose order divides k.
            # A primitive char for p has order m > 1 and m | (p-1).
            count = 0
            for m in get_divisors(k):
                if m == 1:
                    continue
                # Order of primitive character must divide phi(p) = p-1
                if (p - 1) % m == 0:
                    count += phi(m)
            
            num_choices_per_prime.append(count)
            total_prod *= count
        terms[k] = total_prod
        
        # We print the intermediate results for clarity, only for the relevant terms.
        if k in [1, 2, 3, 6]:
            print(f"For each prime in {primes_of_d}, the number of primitive characters with order dividing {k} is:")
            print(f"[{', '.join(map(str, num_choices_per_prime))}]")
            equation_str = " * ".join(map(str, num_choices_per_prime))
            print(f"Total combinations where all orders divide {k}: {equation_str} = {total_prod}\n")


    # Apply inclusion-exclusion principle based on prime factors of g.
    # For g=6=2*3, result = T_6 - T_3 - T_2 + T_1
    result = terms.get(6, 0) - terms.get(3, 0) - terms.get(2, 0) + terms.get(1, 0)
    
    term6 = terms.get(6, 0)
    term3 = terms.get(3, 0)
    term2 = terms.get(2, 0)
    term1 = terms.get(1, 0)

    # Note: Since primitive characters must have order > 1, the count for k=1 is always 0.
    print("Applying the Principle of Inclusion-Exclusion for order 6:")
    final_equation = f"{term6} - {term3} - {term2} + {term1}"
    print(f"The number of characters is {final_equation} = {result}")

    # Return the final answer in the required format
    print(f"\nFinal Answer: {result}")
    
solve()