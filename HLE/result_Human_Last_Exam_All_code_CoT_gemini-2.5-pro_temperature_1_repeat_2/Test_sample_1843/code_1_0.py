import math

def prime_factors(n):
    """Returns a dictionary of prime factors of n."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def get_divisors(n):
    """Returns a set of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return divs

_phi_cache = {}
def phi(n):
    """Computes Euler's totient function."""
    if n in _phi_cache:
        return _phi_cache[n]
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
    _phi_cache[n] = result
    return result

_prim_orders_cache = {}
def get_primitive_orders(p, k):
    """Returns a dictionary {order: count} for primitive characters of conductor p^k."""
    if (p, k) in _prim_orders_cache:
        return _prim_orders_cache[(p, k)]

    orders = {}
    if p == 2:
        if k == 2:  # Conductor 4
            orders = {2: 1}
        elif k == 3:  # Conductor 8
            orders = {2: 2}
        elif k >= 4: # Conductor 2^k
             orders = {2: 1, 2**(k-2): 1}
    else:  # p is an odd prime
        if k == 1:
            group_order = p - 1
            for o in get_divisors(group_order):
                if o > 1:
                    orders[o] = phi(o)
        else:
            group_order = phi(p**k)
            subgroup_order = phi(p**(k-1))
            divs_group = get_divisors(group_order)
            divs_subgroup = get_divisors(subgroup_order)
            for o in divs_group:
                if o not in divs_subgroup:
                    orders[o] = phi(o)
    
    _prim_orders_cache[(p, k)] = orders
    return orders

def count_prim_chars_order_divides_k(p, k, d):
    """Counts primitive characters of conductor p^k whose order divides d."""
    prim_orders = get_primitive_orders(p, k)
    count = 0
    for order, num in prim_orders.items():
        if d % order == 0:
            count += num
    return count

def solve():
    """Main function to solve the problem."""
    N = 36036
    d = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order d = {d}.")
    
    N_factors = prime_factors(N)
    N_factors_str = " * ".join([f"{p}^{k}" for p, k in N_factors.items()])
    print(f"\nStep 1: The prime factorization of N is {N_factors_str}.")

    print("\nStep 2: We use inclusion-exclusion. The number is N(6) - N(3) - N(2) + N(1),")
    print("where N(k) is the number of primitive characters with conductor N and order dividing k.")

    d_prime_factors = list(prime_factors(d).keys())
    
    final_result = 0
    
    terms = {}

    # Iterate through divisors of d for the inclusion-exclusion formula
    divs_of_d = sorted(list(get_divisors(d)), reverse=True)
    
    for k in divs_of_d:
        print(f"\nCalculating N({k}):")
        
        term_count = 1
        term_factors_str = []
        for p, a in N_factors.items():
            count_for_factor = count_prim_chars_order_divides_k(p, a, k)
            print(f"  - For conductor {p**a}, number of primitive characters with order dividing {k} is {count_for_factor}.")
            term_factors_str.append(str(count_for_factor))
            term_count *= count_for_factor
        
        print(f"N({k}) = {' * '.join(term_factors_str)} = {term_count}")
        terms[k] = term_count

    n6 = terms.get(6, 0)
    n3 = terms.get(3, 0)
    n2 = terms.get(2, 0)
    n1 = terms.get(1, 0)
    
    final_result = n6 - n3 - n2 + n1

    print("\nStep 3: Final Calculation")
    print(f"Result = N(6) - N(3) - N(2) + N(1)")
    print(f"Result = {n6} - {n3} - {n2} + {n1} = {final_result}")
    
    print(f"\nThe number of primitive Dirichlet characters of conductor {N} and order {d} is {final_result}.")

solve()