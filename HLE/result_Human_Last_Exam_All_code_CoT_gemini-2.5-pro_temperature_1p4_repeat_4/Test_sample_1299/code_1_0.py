import math

def get_prime_factorization_exponents(n):
    """
    Calculates the exponents of the prime factorization of n.
    For example, if n = 12 = 2^2 * 3^1, it returns [2, 1].
    """
    if n <= 0:
        return []
    if n == 1:
        return []
    
    exponents = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            count = 0
            while temp_n % d == 0:
                count += 1
                temp_n //= d
            exponents.append(count)
        d += 1
    if temp_n > 1:
        exponents.append(1)
    return exponents

def calculate_cardinalities(l):
    """
    Calculates |U_l| and |T_l| for a given integer l.
    """
    print(f"For l = {l}:")
    
    exponents = get_prime_factorization_exponents(l)
    s = len(exponents)

    # --- Part A: Calculate |U_l| ---
    # Formula: 2^max(0, s-1)
    if s == 0: # Corresponds to l=1
        val_A = 1
        print(f"A) |U_{l}| = 2^max(0, {s}-1) = 2^0 = {val_A}")
    else:
        val_A = 2**(s - 1)
        print(f"A) |U_{l}| = 2^({s}-1) = {val_A}")

    # --- Part B: Calculate |T_l| ---
    # Formula: (product of (1+2e_i)) - 1 + 0^(sum of e_i)
    prod = 1
    prod_str = []
    for e in exponents:
        term = 1 + 2 * e
        prod *= term
        prod_str.append(f"(1+2*{e})")
        
    sum_e = sum(exponents)

    val_B = prod - 1 + (0 if sum_e > 0 else 1)

    print("B) |T_l| is calculated from the exponents e_i:", exponents)
    if not exponents: # Case l=1
        print(f"   For l=1, s=0, the product term is 1 and the sum of exponents is 0.")
        print(f"   |T_{l}| = 1 - 1 + 0^0 = 1 - 1 + 1 = {val_B}")
    else:
        print(f"   Product term = {'*'.join(prod_str)} = {prod}")
        print(f"   Sum of exponents = {sum_e}")
        print(f"   |T_{l}| = {prod} - 1 + 0^{sum_e} = {prod} - 1 + 0 = {val_B}")
        
    print("-" * 20)

# Demonstrate with a few examples
calculate_cardinalities(12) # l = 2^2 * 3^1
calculate_cardinalities(1)  # special case l=1
calculate_cardinalities(7)  # prime l
calculate_cardinalities(30) # l = 2*3*5
