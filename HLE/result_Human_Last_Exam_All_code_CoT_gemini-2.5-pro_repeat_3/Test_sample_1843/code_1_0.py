import math

def phi(n):
    """Computes Euler's totient function."""
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

def get_prime_power(m):
    """Decomposes m into p^a."""
    if m == 1:
        return 1, 0
    num = m
    p = 2
    while p * p <= num:
        if num % p == 0:
            a = 0
            while num % p == 0:
                num //= p
                a += 1
            if num == 1:
                return p, a
        p += 1
    if num > 1:
        return num, 1
    return None, None # Should not happen for p^a

def count_primitive_chars_order_k(m, k):
    """
    Computes the number of primitive Dirichlet characters of conductor m and order k.
    """
    if k == 1:
        # Primitive characters cannot have order 1 (only the principal character has order 1, conductor 1)
        return 0

    p, a = get_prime_power(m)

    if p == 2:
        if a == 1: # m = 2
            return 0
        if a == 2: # m = 4
            return 1 if k == 2 else 0
        if a >= 3:
            # Not needed for N=36036, but for completeness:
            # G_2^a is C_2 x C_{2^{a-2}}. G_{2^{a-1}} is C_2 x C_{2^{a-3}}
            # This part is more complex and not implemented as it is not required.
            return 0
    else: # p is odd
        if a == 1: # m = p
            if k > 1 and (p - 1) % k == 0:
                return phi(k)
            else:
                return 0
        else: # m = p^a, a > 1
            phi_pa = (p - 1) * (p**(a - 1))
            phi_pa_prev = (p - 1) * (p**(a - 2))
            
            count_pa = phi(k) if phi_pa % k == 0 else 0
            count_pa_prev = phi(k) if phi_pa_prev % k == 0 else 0
            return count_pa - count_pa_prev

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of conductor N=36036 and order 6.
    """
    N = 36036
    order_k = 6
    
    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order {order_k}.")
    print("N factors as 2^2 * 3^2 * 7^1 * 11^1 * 13^1 = 4 * 9 * 7 * 11 * 13.")
    print("We use the principle of inclusion-exclusion.\n")

    moduli = [4, 9, 7, 11, 13]
    
    divisors_of_6 = [1, 2, 3, 6]
    divisors_of_3 = [1, 3]
    divisors_of_2 = [1, 2]

    counts = {
        'div_6': [],
        'div_3': [],
        'div_2': [],
        'div_1': []
    }

    for m in moduli:
        count_div_6 = sum(count_primitive_chars_order_k(m, k) for k in divisors_of_6)
        count_div_3 = sum(count_primitive_chars_order_k(m, k) for k in divisors_of_3)
        count_div_2 = sum(count_primitive_chars_order_k(m, k) for k in divisors_of_2)
        count_div_1 = count_primitive_chars_order_k(m, 1)
        
        counts['div_6'].append(count_div_6)
        counts['div_3'].append(count_div_3)
        counts['div_2'].append(count_div_2)
        counts['div_1'].append(count_div_1)
    
    total_div_6 = math.prod(counts['div_6'])
    total_div_3 = math.prod(counts['div_3'])
    total_div_2 = math.prod(counts['div_2'])
    total_div_1 = math.prod(counts['div_1']) # This will be 0

    print("Number of primitive characters whose order divides 6 for conductors 4, 9, 7, 11, 13:")
    print(f"Counts: {counts['div_6']}")
    print(f"Product = {' * '.join(map(str, counts['div_6']))} = {total_div_6}\n")
    
    print("Number of primitive characters whose order divides 3 for conductors 4, 9, 7, 11, 13:")
    print(f"Counts: {counts['div_3']}")
    print(f"Product = {' * '.join(map(str, counts['div_3']))} = {total_div_3}\n")

    print("Number of primitive characters whose order divides 2 for conductors 4, 9, 7, 11, 13:")
    print(f"Counts: {counts['div_2']}")
    print(f"Product = {' * '.join(map(str, counts['div_2']))} = {total_div_2}\n")

    result = total_div_6 - total_div_3 - total_div_2 + total_div_1
    
    print("Final calculation using inclusion-exclusion:")
    print(f"Result = (chars order | 6) - (chars order | 3) - (chars order | 2) + (chars order | 1)")
    print(f"Result = {total_div_6} - {total_div_3} - {total_div_2} + {total_div_1} = {result}")

    print("\nEach number in the final equation:")
    # Using 'x' for multiplication symbol for clarity
    term1 = ' x '.join(map(str, counts['div_6']))
    term2 = ' x '.join(map(str, counts['div_3']))
    term3 = ' x '.join(map(str, counts['div_2']))
    print(f"({term1}) - ({term2}) - ({term3}) + {total_div_1} = {result}")

solve()
<<<100>>>