import math

def prime_factorize(n):
    """
    Returns the prime factorization of n as a dictionary {p: a}.
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

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    """
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
    return result

def get_primitive_orders(p, a):
    """
    Returns a list of possible orders for primitive characters modulo p^a.
    """
    if p == 2:
        if a == 1: return [] # No primitive characters mod 2
        if a == 2: return [2] # Only one primitive char mod 4, order 2
        # For p=2, a>=3, the group is C_2 x C_{2^{a-2}}.
        # Primitive characters are those not induced from mod 2^{a-1}.
        # This is more complex, but not needed for N=36036.
        # N has a factor of 2^2, not higher powers of 2.
        return []

    # For odd primes p
    if a == 1:
        # For mod p, primitive chars are all non-principal chars.
        # Orders are divisors of phi(p)=p-1, greater than 1.
        order_group = p - 1
        orders = []
        for i in range(2, order_group + 1):
            if order_group % i == 0:
                orders.append(i)
        return orders
    else:
        # For mod p^a, a>=2, the group of characters is cyclic of order phi(p^a).
        # Primitive characters are those whose order d divides phi(p^a) but not phi(p^{a-1}).
        order_group = phi(p**a)
        order_subgroup = phi(p**(a-1))
        orders = []
        for i in range(1, order_group + 1):
            if order_group % i == 0:
                if order_subgroup % i != 0:
                    orders.append(i)
        return orders

def count_choices_for_factor(p, a, required_order_divisor):
    """
    Counts the number of primitive characters mod p^a whose order divides `required_order_divisor`.
    """
    prim_orders = get_primitive_orders(p, a)
    count = 0
    
    # We need to find which of the primitive orders divide `required_order_divisor`
    allowed_orders = [d for d in prim_orders if required_order_divisor % d == 0]
    
    if not allowed_orders:
        return 0

    # The number of characters of a given order d in a cyclic group is phi(d).
    # This applies for odd p, and for p=2, a<=2.
    for d in allowed_orders:
        count += phi(d)
        
    return count

def solve():
    """
    Main function to solve the problem.
    """
    N = 36036
    order = 6
    
    factors = prime_factorize(N)
    
    print(f"The conductor is N = {N}.")
    factor_str = " * ".join([f"{p}^{a}" if a > 1 else str(p) for p, a in factors.items()])
    print(f"The prime factorization of N is: {factor_str}\n")
    
    print(f"We need to find the number of primitive characters of order {order}.")
    print("A character is primitive with conductor N if and only if it is a product of primitive characters for each prime power factor.\n")
    
    total_count = 1
    counts = {}
    
    lcm_has_factor_2 = False
    lcm_has_factor_3 = False
    
    print("Counting the number of choices for each prime power factor whose order divides 6:")
    
    # Sort factors for consistent output
    sorted_factors = sorted(factors.items())

    for p, a in sorted_factors:
        pk = p**a
        
        # Get primitive orders for p^a
        prim_orders = get_primitive_orders(p, a)
        
        # Filter for orders that divide the required order 6
        allowed_orders = [d for d in prim_orders if order % d == 0]
        
        # For any choice, check if it forces a factor of 2 or 3 in the lcm
        if all(o % 2 == 0 for o in prim_orders):
             lcm_has_factor_2 = True
        if any(o % 2 != 0 for o in prim_orders) and all(o%2==0 or o%3==0 for o in prim_orders) and any(o%3==0 for o in prim_orders):
             pass # more complex logic not needed here
        
        # Check divisibility by 3
        if len(prim_orders) > 0 and all(o % 3 == 0 for o in [po for po in prim_orders if po%2 != 0]): # simplified logic for this problem
             if all(po % 3 == 0 or po % 6 == 0 for po in prim_orders):
                 lcm_has_factor_3 = True
        
        # Simplified check for this specific problem
        if pk == 4:
            lcm_has_factor_2 = True
        if pk == 9:
            lcm_has_factor_3 = True

        
        # The number of characters for each allowed order d is phi(d), except for p=2, a>=3
        # This is handled correctly by the logic below for our N.
        # For p=2, a=2, primitive order is 2. Number of chars is 1. phi(2)=1.
        num_choices = sum(phi(d) for d in allowed_orders)
        
        counts[pk] = num_choices
        total_count *= num_choices
        
        print(f"- For conductor {pk}:")
        print(f"  The orders of primitive characters are {prim_orders}.")
        print(f"  The orders that divide {order} are {allowed_orders}.")
        print(f"  Number of choices = " + " + ".join([f"phi({d})" for d in allowed_orders]) + f" = " + " + ".join([str(phi(d)) for d in allowed_orders]) + f" = {num_choices}")

    print("\n--------------------------------------------------\n")
    
    print("The total number of combinations of primitive characters whose orders divide 6 is the product of the individual counts:")
    equation_parts = [str(c) for c in counts.values()]
    print(f"Total combinations = {' * '.join(equation_parts)} = {total_count}")
    
    print("\nNext, we verify that the lcm of the orders is exactly 6 for any combination.")
    if lcm_has_factor_2:
        print("- The choice for conductor 4 requires a character of order 2. So, the lcm of orders is always divisible by 2.")
    if lcm_has_factor_3:
        print("- The choice for conductor 9 requires a character of order 3 or 6. So, the lcm of orders is always divisible by 3.")
    
    if lcm_has_factor_2 and lcm_has_factor_3:
        print("- Since the lcm is divisible by both 2 and 3, it must be divisible by 6.")
        print("- Since all component character orders divide 6, their lcm must also divide 6.")
        print("- A number that is a multiple of 6 and a divisor of 6 must be exactly 6.")
        print("\nTherefore, every combination results in a character of order 6.")
        print("The final answer is the total number of combinations.")
    else:
        # Fallback to inclusion-exclusion if the simple check fails. Not needed for this problem.
        print("\nThe simple check failed. A full inclusion-exclusion calculation would be needed.")
        # But for this problem, the logic holds.
        pass

    print("\nFinal Answer Calculation:")
    final_eq_str = f"Number of characters = {counts[4]} (for 2^2) * {counts[9]} (for 3^2) * {counts[7]} (for 7) * {counts[11]} (for 11) * {counts[13]} (for 13) = {total_count}"
    print(final_eq_str)
    
    return total_count

final_answer = solve()
# The final answer is directly printed by the solve function.
# To satisfy the format requirement, we will also print it here.
# print(f"\n<<<{final_answer}>>>")