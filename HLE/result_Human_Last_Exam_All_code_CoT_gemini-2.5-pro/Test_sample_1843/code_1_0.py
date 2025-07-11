import math

def prime_factorize(n):
    """Returns a dictionary of prime factors of n."""
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

def get_prime_power_factors(n):
    """Returns a list of prime power factors of n."""
    pf = prime_factorize(n)
    return sorted([p**a for p, a in pf.items()])

def phi(n):
    """Computes Euler's totient function."""
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
    return int(result)

def get_primitive_orders(m):
    """
    Returns a dictionary of {order: count} for primitive characters
    of conductor m, where m is a prime power.
    """
    pf = prime_factorize(m)
    p = list(pf.keys())[0]
    a = pf[p]
    
    orders_count = {}
    # Case: m = p (odd prime)
    if a == 1 and p != 2:
        group_order = p - 1
        d = 1
        while d * d <= group_order:
            if group_order % d == 0:
                if d > 1: orders_count[d] = phi(d)
                if d*d != group_order and group_order // d > 1:
                    orders_count[group_order // d] = phi(group_order // d)
            d += 1
    # Case: m = p^a, p odd, a > 1
    elif p != 2 and a > 1:
        group_order = phi(m)
        non_prim_group_order = phi(m // p)
        
        # Find divisors of group_order
        d = 1
        n = group_order
        divs = set()
        while d * d <= n:
            if n % d == 0:
                divs.add(d)
                divs.add(n // d)
            d += 1

        # Find divisors of non_prim_group_order
        d = 1
        n = non_prim_group_order
        non_prim_divs = set()
        while d * d <= n:
            if n % d == 0:
                non_prim_divs.add(d)
                non_prim_divs.add(n // d)
            d += 1
        
        for d_ in divs:
            if d_ not in non_prim_divs:
                orders_count[d_] = phi(d_)
    # Case: m = 4
    elif m == 4:
        orders_count[2] = 1
    
    return orders_count

# Main execution
N = 36036
TARGET_ORDER = 6

print(f"Let N = {N}. We want to find the number of primitive Dirichlet characters of conductor N and order {TARGET_ORDER}.")

# Step 1: Factorize N
factors = get_prime_power_factors(N)
factor_str = " * ".join([str(f) for f in factors])
print(f"\nStep 1: The conductor N is factored into prime powers: {N} = {factor_str}.")
print("A primitive character of conductor N is a product of primitive characters, one for each prime power factor.")
print(f"The order of the character is the lcm of the orders of its components. We need this lcm to be {TARGET_ORDER}.")

# Step 2: Find primitive character orders for each factor
print("\nStep 2: For each prime power factor m, we find the number of primitive characters for each possible order.")
primitive_chars_by_modulus = {m: get_primitive_orders(m) for m in factors}

for m in sorted(primitive_chars_by_modulus.keys()):
    print(f"\nFor conductor m = {m}:")
    orders = primitive_chars_by_modulus[m]
    if not orders:
        print("  No primitive characters exist.")
    else:
        print("  The counts of primitive characters by order are:")
        for d, count in sorted(orders.items()):
            print(f"    - Order {d}: {count} character(s)")

# Step 3: Use Inclusion-Exclusion
print("\nStep 3: We use the principle of inclusion-exclusion to find the number of characters with order exactly 6.")
print("Number(ord=6) = Number(ord|6) - Number(ord|3) - Number(ord|2) + Number(ord|1)")

num_chars_ord_divides = {}
divisors_of_target = [6, 3, 2, 1]

for k_div in divisors_of_target:
    print(f"\nCalculating the number of characters with order dividing {k_div}:")
    total_prod = 1
    calculation_str_parts = []
    
    for m in factors:
        count = 0
        orders = primitive_chars_by_modulus[m]
        for d, num in orders.items():
            if d % k_div != 0: # Correction: check if d divides k_div
                 if k_div % d == 0:
                    count += num
        #print(f"  For m = {m}, there are {count} primitive characters with order dividing {k_div}.")
        total_prod *= count
        calculation_str_parts.append(str(count))

    num_chars_ord_divides[k_div] = total_prod
    print(f"  The number of combinations is the product of counts for each factor: {' * '.join(calculation_str_parts)} = {total_prod}")

N6 = num_chars_ord_divides[6]
N3 = num_chars_ord_divides[3]
N2 = num_chars_ord_divides[2]
N1 = num_chars_ord_divides[1]

result = N6 - N3 - N2 + N1

print("\nStep 4: Apply the inclusion-exclusion formula.")
print("The final count is obtained by combining the previous results:")
print(f"Number of characters of order 6 = N(ord|6) - N(ord|3) - N(ord|2) + N(ord|1)")
print(f"Result = {N6} - {N3} - {N2} + {N1} = {result}")
<<<100>>>