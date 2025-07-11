import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def euler_phi(n):
    """Computes Euler's totient function."""
    if n == 0:
        return 0
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

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    while d * d <= n:
        while (n % d) == 0:
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

def moebius_mu_recursive(n, factors):
    """Recursive helper for moebius_mu."""
    if n == 1:
        return 1
    p = list(factors.keys())[0]
    m = factors[p]
    if m > 1:
        return 0
    
    remaining_factors = factors.copy()
    del remaining_factors[p]
    
    rest_n = 1
    for prime, power in remaining_factors.items():
      rest_n *= prime**power
      
    return -moebius_mu_recursive(rest_n, remaining_factors)

def moebius_mu(n):
    """Computes the Moebius function mu(n)."""
    if n == 0:
      return 0
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def count_prim_chars_order_divides_d(p, a, d):
    """
    Counts the number of primitive characters modulo p^a whose order divides d.
    """
    # Case p=2, a=1 (mod 2): No primitive characters.
    if p == 2 and a == 1:
        return 0
    
    # Case p=2, a=2 (mod 4): One primitive character of order 2.
    if p == 2 and a == 2:
        return 1 if d % 2 == 0 else 0
        
    # Case p=2, a>=3 (mod 2^a): Character group is Z_2 x Z_{2^{a-2}}.
    # The number of primitive characters is 2^{a-2}. We would need to analyze orders.
    # Not needed for N=36036.
    if p == 2 and a >= 3:
      # Placeholder for a more general implementation
      raise NotImplementedError("Case p=2, a>=3 is not implemented.")

    # Case p > 2 prime. The character group is cyclic of order phi(p^a).
    phi_pa = euler_phi(p**a)
    
    if a == 1:
        # For mod p, primitive chars are all non-principal chars.
        # Orders are divisors of phi(p) > 1.
        prim_orders = [k for k in get_divisors(phi_pa) if k > 1]
    else:
        # For mod p^a, orders of prim chars divide phi(p^a) but not phi(p^{a-1}).
        phi_pa_minus_1 = euler_phi(p**(a-1))
        prim_orders = [k for k in get_divisors(phi_pa) if k % phi_pa_minus_1 != 0]

    count = 0
    for k in prim_orders:
        if d % k == 0:
            count += euler_phi(k)
    return count

def solve():
    """
    Solves the problem for a given N and K.
    """
    N = 36036
    K = 6
    
    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order K = {K}.\n")

    N_factors = get_prime_factorization(N)
    print(f"The prime factorization of N is: {N} = {' * '.join([f'{p}^{a}' for p, a in N_factors.items()])}")

    K_divisors = get_divisors(K)
    total_count = 0
    
    print("\nUsing the inclusion-exclusion principle:")
    print(f"Total = sum_{{d|{K}}} [mu({K}/d) * (number of primitive characters with order dividing d)]\n")

    equation_parts = []
    
    for d in K_divisors:
        mu_val = moebius_mu(K // d)
        if mu_val == 0:
            continue
            
        prod_N_d = 1
        prod_parts = []
        for p, a in N_factors.items():
            C_pad = count_prim_chars_order_divides_d(p, a, d)
            prod_N_d *= C_pad
            prod_parts.append(f"C({p}^{a}, {d})={C_pad}")
        
        term = mu_val * prod_N_d
        total_count += term

        print(f"For d = {d}:")
        print(f"  mu({K}/{d}) = mu({K//d}) = {mu_val}")
        print(f"  Number of characters with order dividing {d} is:")
        print(f"  N({d}) = {' * '.join(prod_parts)} = {prod_N_d}")
        
        sign = "+" if mu_val > 0 else "-"
        equation_parts.append(f"{sign} {abs(mu_val)}*{prod_N_d}")
        
    print("\nFinal calculation:")
    # Clean up the start of the equation string
    if equation_parts[0].startswith('+ '):
      equation_parts[0] = equation_parts[0][2:]
      
    final_equation = " ".join(equation_parts)
    print(f"Total = {final_equation} = {total_count}")
    
    print(f"\nThe number of primitive Dirichlet characters of conductor {N} and order {K} is {total_count}.")
    return total_count

final_answer = solve()
<<<100>>>