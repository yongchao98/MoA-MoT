import math
from collections import defaultdict

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return dict(factors)

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

memo_phi = {}
def phi(n):
    """Euler's totient function."""
    if n in memo_phi:
        return memo_phi[n]
    if n == 1:
        memo_phi[1] = 1
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
    memo_phi[n] = result
    return result

memo_mobius = {}
def mobius(n):
    """Mobius function mu(n)."""
    if n in memo_mobius:
        return memo_mobius[n]
    if n == 1:
        memo_mobius[1] = 1
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            memo_mobius[n] = 0
            return 0
    result = (-1)**len(factors)
    memo_mobius[n] = result
    return result

def count_primitive_orders(p, a):
    """
    Counts the number of primitive Dirichlet characters modulo p^a by their order.
    Returns a dictionary {order: count}.
    """
    if p == 2:
        if a == 1: return {}
        if a == 2: return {2: 1}
        # For a >= 3, all 2^(a-2) primitive characters have order 2^(a-2).
        return {2**(a - 2): 2**(a - 2)}
    
    # p is an odd prime
    if a == 1:
        m = p - 1
        orders = {}
        divs = get_divisors(m)
        for k in divs:
            if k > 1:
                orders[k] = phi(k)
        return orders
    else: # a >= 2
        m = phi(p**a)
        m_non_primitive = phi(p**(a - 1))
        
        orders = {}
        divs_m = get_divisors(m)
        for k in divs_m:
            if k > 0 and m_non_primitive % k != 0:
                orders[k] = phi(k)
        return orders

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of conductor N and order K.
    """
    N = 36036
    K = 6

    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order K = {K}.")
    
    # Step 1: Prime factorization of N
    factors_N = get_prime_factorization(N)
    print(f"\nStep 1: The prime factorization of N is {N} = {' * '.join([f'{p}^{a}' for p, a in factors_N.items()])}.")

    # Step 2: Get orders for primitive characters for each factor
    print("\nStep 2: For each prime power factor p^a, we count the primitive characters by their order.")
    char_orders_by_factor = {}
    for p, a in sorted(factors_N.items()):
        pa_val = p**a
        char_orders_by_factor[pa_val] = count_primitive_orders(p, a)
        print(f"  - For conductor {pa_val}: {char_orders_by_factor[pa_val]}")

    # Step 3: Use Mobius inversion. Calculate N(d) for d|K.
    divs_K = get_divisors(K)
    N_d = {}

    print(f"\nStep 3: We use Mobius inversion. We need to compute N(d) for each divisor d of K={K}.")
    print(f"The divisors of {K} are {divs_K}.")

    for d in divs_K:
        num_product = 1
        print(f"\n  Calculating N({d}), the number of characters with order dividing {d}:")
        for pa, orders in sorted(char_orders_by_factor.items()):
            count_d = 0
            for order, num in orders.items():
                if d % order == 0:
                    count_d += num
            print(f"    - For conductor {pa}, there are {count_d} primitive characters with order dividing {d}.")
            num_product *= count_d
        N_d[d] = num_product
        print(f"  N({d}) = The product of these counts = {num_product}")

    # Step 4: Apply the Mobius inversion formula
    print(f"\nStep 4: Apply the Mobius inversion formula: Result = sum(mu(K/d) * N(d)) for d|K.")
    total_count = 0
    equation_parts = []
    
    # Building the equation string
    for d in divs_K:
        mu_val = mobius(K // d)
        n_val = N_d[d]
        if n_val == 0:
            continue # Skip terms that are zero
        
        sign = "+"
        if mu_val < 0:
            sign = "-"
        
        if len(equation_parts) > 0:
            equation_parts.append(f" {sign} ")

        equation_parts.append(f"{abs(mu_val)} * {n_val}")

    # Calculate final result
    for d in divs_K:
        mu_val = mobius(K // d)
        total_count += mu_val * N_d[d]

    # Combine equation parts with the actual numbers
    final_equation_str = ""
    for d in divs_K:
      mu_val = mobius(K//d)
      term = mu_val * N_d[d]
      if final_equation_str == "":
        final_equation_str += f"({mu_val}) * {N_d[d]}"
      else:
        final_equation_str += f" + ({mu_val}) * {N_d[d]}"

    print(f"Result = {final_equation_str}")
    print(f"Result = {total_count}")
    
    print(f"\nThe number of primitive Dirichlet characters of conductor {N} and order {K} is {total_count}.")

solve()
<<<100>>>