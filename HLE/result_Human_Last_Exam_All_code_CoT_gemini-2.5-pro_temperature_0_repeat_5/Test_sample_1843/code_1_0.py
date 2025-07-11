import math

def phi(n):
    """Computes Euler's totient function."""
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
    return result

def count_order_in_cyclic(d, m):
    """Counts the number of elements of order d in a cyclic group of order m."""
    if m % d == 0:
        return phi(d)
    return 0

def get_primitive_char_counts_by_order(p, a):
    """
    Calculates the number of primitive characters for conductor p^a, grouped by order.
    Returns a dictionary {order: count}.
    """
    counts = {}
    if p == 2:
        if a == 1: return {}
        if a == 2: return {2: 1} # Primitive char mod 4 has order 2
        # For a >= 3, group is C2 x C_{2^(a-2)}, imprimitive is C2 x C_{2^(a-3)}
        # This is more complex, but not needed for N=36036
        return {} 
    else: # p is odd
        group_order = phi(p**a)
        if a == 1:
            # All non-trivial characters are primitive
            for d in range(2, group_order + 1):
                if group_order % d == 0:
                    counts[d] = phi(d)
        else: # a > 1
            imprimitive_group_order = phi(p**(a-1))
            # Iterate through all possible orders d in the full group
            for d in range(2, group_order + 1):
                if group_order % d == 0:
                    num_total = phi(d)
                    num_imprimitive = count_order_in_cyclic(d, imprimitive_group_order)
                    num_primitive = num_total - num_imprimitive
                    if num_primitive > 0:
                        counts[d] = num_primitive
    return counts

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of a given conductor and order.
    """
    N = 36036
    target_order = 6
    
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order {target_order}.")
    print("A primitive character chi mod N can be uniquely decomposed into a product of primitive characters chi_{p^a} for each prime power p^a in the factorization of N.")
    print("The order of chi is the least common multiple (lcm) of the orders of its components chi_{p^a}.")
    
    # Step 1: Factor N
    factors = {2: 2, 3: 2, 7: 1, 11: 1, 13: 1} # N = 4 * 9 * 7 * 11 * 13
    print(f"\nStep 1: The prime power factorization of N is N = 2^2 * 3^2 * 7^1 * 11^1 * 13^1.")

    # Step 2: For each component, find the number of primitive characters whose order divides the target order.
    print(f"\nStep 2: For a character to have order {target_order}, the order of each of its components must divide {target_order}.")
    print("We count the number of primitive characters for each prime power factor whose order divides 6.")

    component_counts = {}
    
    # Analysis for each factor
    for p, a in sorted(factors.items()):
        pa = p**a
        prim_chars = get_primitive_char_counts_by_order(p, a)
        count_pa = sum(c for d, c in prim_chars.items() if target_order % d == 0)
        component_counts[pa] = count_pa
        
        valid_orders = [d for d in prim_chars if target_order % d == 0]
        print(f"For conductor {pa}, primitive characters with orders dividing 6 are: {valid_orders}. This gives {count_pa} choice(s).")

    # Step 3: Check if the lcm is always 6
    print("\nStep 3: Check if the lcm of the orders of these components is always 6.")
    print("The order of the component for conductor 4 is 2, so the final order is a multiple of 2.")
    print("The order of the component for conductor 9 is 3 or 6, so the final order is a multiple of 3.")
    print("Thus, the final order is a multiple of lcm(2, 3) = 6.")
    print("Since all component orders divide 6, their lcm must also divide 6.")
    print("Therefore, any combination of these selected components results in a character of order exactly 6.")

    # Step 4: Calculate the total number
    print("\nStep 4: The total number of such characters is the product of the number of choices for each component.")
    
    total_count = 1
    equation_parts = []
    for pa in sorted(component_counts.keys()):
        count = component_counts[pa]
        total_count *= count
        equation_parts.append(str(count))
        
    print(f"Number of primitive characters of order 6 = {' * '.join(equation_parts)} = {total_count}")

solve()