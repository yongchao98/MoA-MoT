import math

def prime_factorize(n):
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
    """Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return int(result)

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_character_count():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    N = 36036
    ORDER = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order {ORDER}.")
    print("-" * 70)

    # Step 1: Factorize N
    print("Step 1: Factorize the conductor N.")
    factors = prime_factorize(N)
    factor_str = " * ".join([f"{p}^{a}" if a > 1 else str(p) for p, a in sorted(factors.items())])
    print(f"N = {N} = {factor_str}")
    print("\nA character chi mod N is primitive with conductor N iff it is a product of primitive")
    print("characters chi_{p^a} with conductor p^a for each prime power factor p^a of N.")
    print(f"The order of chi must be lcm(orders of components) = {ORDER}.")
    print("-" * 70)

    print(f"Step 2: For each factor p^a, count primitive characters with order dividing {ORDER}.")
    
    total_count = 1
    component_counts = {}
    
    # Analysis for each prime power factor
    # Conductor 4 = 2^2
    p, a = 2, 2
    pa = p**a
    print(f"\nConductor {pa} ({p}^{a}):")
    group_order = phi(pa) if a > 1 else phi(p)
    print(f"  Group of characters is C_{group_order}. Primitive characters have order 2.")
    count_4 = phi(2)
    print(f"  Orders dividing {ORDER}: [2]. Count = phi(2) = {count_4}.")
    component_counts[pa] = count_4
    
    # Conductor 9 = 3^2
    p, a = 3, 2
    pa = p**a
    print(f"\nConductor {pa} ({p}^{a}):")
    group_order = phi(pa)
    print(f"  Group of characters is C_{group_order}. Primitive characters have orders 3, 6.")
    count_9 = phi(3) + phi(6)
    print(f"  Orders dividing {ORDER}: [3, 6]. Count = phi(3) + phi(6) = {phi(3)} + {phi(6)} = {count_9}.")
    component_counts[pa] = count_9

    # Prime conductors
    for p in [7, 11, 13]:
        pa = p
        print(f"\nConductor {pa}:")
        group_order = phi(p)
        prim_orders = [d for d in get_divisors(group_order) if d > 1]
        print(f"  Group of characters is C_{group_order}. Primitive orders: {prim_orders}.")
        
        orders_dividing_k = [d for d in prim_orders if ORDER % d == 0]
        count = sum(phi(d) for d in orders_dividing_k)
        
        phi_calcs = " + ".join([f"phi({d})" for d in orders_dividing_k])
        phi_vals = " + ".join([str(phi(d)) for d in orders_dividing_k])
        print(f"  Orders dividing {ORDER}: {orders_dividing_k}. Count = {phi_calcs} = {phi_vals} = {count}.")
        component_counts[pa] = count

    print("-" * 70)
    print("Step 3: Combine results and verify order condition.")
    
    # Sort component_counts by key for consistent output
    sorted_counts = {k: v for k, v in sorted(component_counts.items())}
    
    calc_str = " * ".join(map(str, sorted_counts.values()))
    total_count = math.prod(sorted_counts.values())
    print(f"Total combinations of component characters = {calc_str} = {total_count}")

    print("\nVerification of LCM:")
    print("  - Character for conductor 4 has order 2.")
    print("  - Character for conductor 9 has order 3 or 6.")
    print("  - The lcm is guaranteed to be divisible by lcm(2, 3) = 6.")
    print("  - All component orders divide 6, so their lcm must divide 6.")
    print("  => The final order is always exactly 6.")
    print("-" * 70)

    print("The number of primitive characters of conductor 36036 and order 6 is:")
    final_eq = f"{sorted_counts[4]} * {sorted_counts[9]} * {sorted_counts[7]} * {sorted_counts[11]} * {sorted_counts[13]} = {total_count}"
    print(final_eq)
    
    return total_count

final_answer = solve_character_count()
<<<100>>>