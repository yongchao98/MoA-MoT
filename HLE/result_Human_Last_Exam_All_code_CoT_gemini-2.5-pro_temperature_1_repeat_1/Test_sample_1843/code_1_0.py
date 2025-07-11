import math

def phi(n):
    """Euler's totient function."""
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

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def get_prime_power(n):
    """Finds p, k for n=p^k."""
    if n == 4: return 2, 2
    temp_n = n
    p = 2
    if temp_n % p == 0:
        k = 0
        while temp_n % p == 0:
            temp_n //= p
            k += 1
        if temp_n == 1: return p, k
    p = 3
    while p * p <= n:
        if n % p == 0:
            k = 0
            temp_n = n
            while temp_n % p == 0:
                temp_n //= p
                k += 1
            if temp_n == 1: return p, k
        p += 2
    return n, 1 # n is prime

def get_primitive_character_orders(modulus):
    """
    Returns a dictionary mapping order to the count of primitive characters
    for a given prime power modulus.
    """
    if modulus == 1:
        return {}
    if modulus == 4:
        # The character group mod 4 is Z_2. The non-trivial character is primitive and has order 2.
        return {2: 1}
    
    p, k = get_prime_power(modulus)
    
    # For p^k, p odd or p=2, k>2, the char group is isomorphic to (Z/p^k Z)*
    # For p odd, it's cyclic of order phi(p^k).
    # For a prime p, primitive characters are all non-principal characters.
    if k == 1:
        group_order = phi(p) # p-1
        orders = {}
        for d in get_divisors(group_order):
            if d > 1:
                orders[d] = phi(d)
        return orders

    # For p^k, k>1, p odd, primitive characters have orders d | phi(p^k) but d does not divide phi(p^{k-1}).
    if p > 2 and k > 1:
        group_order = phi(modulus)
        subgroup_order = phi(modulus // p)
        divs_pk = get_divisors(group_order)
        divs_pk_minus_1 = get_divisors(subgroup_order)
        orders = {}
        for d in divs_pk:
            if d not in divs_pk_minus_1:
                orders[d] = phi(d)
        return orders
    return {}


def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of conductor N and order 6.
    """
    N = 36036
    target_order = 6
    
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order {target_order}.")
    print("\nStep 1: Prime factorization of N.")
    print(f"N = {N} = 2^2 * 3^2 * 7 * 11 * 13 = 4 * 9 * 7 * 11 * 13.")
    
    print("\nStep 2: A character is primitive mod N if it's a product of primitive characters for each prime power factor.")
    print(f"The order of the product character must be {target_order}.")
    
    print(f"\nStep 3: For the lcm of orders to be {target_order}, the order of each component character must divide {target_order}.")
    
    print("\nStep 4: Count the number of suitable primitive characters for each prime power factor.")
    
    moduli = [4, 9, 7, 11, 13]
    total_count = 1
    
    counts = {}

    for p_k in moduli:
        print(f"\n--- Analyzing modulus p^k = {p_k} ---")
        primitive_orders = get_primitive_character_orders(p_k)
        print(f"The possible orders for primitive characters modulo {p_k} are: {list(primitive_orders.keys())}")
        
        count_for_pk = 0
        orders_dividing_target = []
        for o, num in primitive_orders.items():
            if target_order % o == 0:
                count_for_pk += num
                orders_dividing_target.append((o, num))
        
        counts[p_k] = {
            'choices': count_for_pk,
            'details': orders_dividing_target
        }
        total_count *= count_for_pk
        
        print(f"The orders that divide {target_order} are: {[o for o, n in orders_dividing_target]}")
        detail_str = " + ".join([f"{n} (of order {o})" for o, n in orders_dividing_target])
        print(f"Number of suitable primitive characters for modulus {p_k} is {detail_str or '0'} = {count_for_pk}.")

    print("\nStep 5: Verify that the lcm of the orders is exactly 6.")
    o4_orders = [o for o, n in counts[4]['details']]
    o9_orders = [o for o, n in counts[9]['details']]
    print(f"From our analysis, the order for mod 4 must be in {o4_orders}, so the lcm is a multiple of 2.")
    print(f"The order for mod 9 must be in {o9_orders}, so the lcm is a multiple of 3.")
    print("Therefore, the lcm must be a multiple of lcm(2, 3) = 6.")
    print("Since we already ensured each component's order divides 6, the final lcm must be exactly 6.")

    print("\nStep 6: Calculate the total number of characters.")
    print("The total number is the product of the number of choices for each modulus.")
    
    counts_list = [counts[p_k]['choices'] for p_k in moduli]
    calculation_str = " * ".join(map(str, counts_list))
    print(f"Number of characters = {calculation_str} = {total_count}")

solve()