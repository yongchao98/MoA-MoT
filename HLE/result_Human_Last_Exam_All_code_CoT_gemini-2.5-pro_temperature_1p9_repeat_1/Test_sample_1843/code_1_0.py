import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    """
    if not hasattr(phi, "memo"):
        phi.memo = {}
    if n in phi.memo:
        return phi.memo[n]
    
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
    phi.memo[n] = result
    return result

def get_divisors(n):
    """
    Returns a sorted list of all divisors of n.
    """
    if not hasattr(get_divisors, "memo"):
        get_divisors.memo = {}
    if n in get_divisors.memo:
        return get_divisors.memo[n]
        
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    get_divisors.memo[n] = sorted(list(divs))
    return get_divisors.memo[n]

def count_prim_chars_order_divides_d(p, a, d):
    """
    Counts the number of primitive characters modulo p^a whose order divides d.
    """
    # For prime modulus p, primitive characters are the non-trivial ones.
    # The group of characters is cyclic of order p-1.
    if a == 1:
        if p == 2: # Mod 2, only trivial char (conductor 1), so no primitive chars.
            return 0
        group_order = p - 1
        common_order_divs = math.gcd(d, group_order)
        divs = get_divisors(common_order_divs)
        count = 0
        for k in divs:
            if k > 1:  # Order > 1 for non-trivial characters
                count += phi(k)
        return count

    # For modulus p^a, a > 1
    if p == 2:
        if a == 2: # Mod 4, one primitive character of order 2.
            return 1 if d % 2 == 0 else 0
        else: # Mod 2^a for a >= 3 not needed for this problem.
            raise ValueError(f"Case p=2, a={a} is not handled.")
            
    # For odd prime p, a >= 2. Primitive characters are those not induced
    # from modulus p^(a-1).
    phi_pa = phi(p**a)
    phi_pa_minus_1 = phi(p**(a-1))

    count = 0
    # Sum phi(k) for all k that divide d, divide phi(p^a), but NOT divide phi(p^(a-1)).
    for k in get_divisors(d):
        if k > 0 and phi_pa % k == 0 and phi_pa_minus_1 % k != 0:
            count += phi(k)
            
    return count

def solve_problem():
    """
    Main function to solve the problem and print the detailed steps.
    """
    N = 36036
    factors = [(2, 2), (3, 2), (7, 1), (11, 1), (13, 1)]
    
    print(f"Let N = {N}. We want to find the number of primitive Dirichlet characters of conductor N and order 6.")
    print("Step 1: Factorize N.")
    print("N = 36036 = 2^2 * 3^2 * 7 * 11 * 13")
    
    print("\nStep 2: Use the formula based on inclusion-exclusion for order 6.")
    print("Num(order 6) = N_6 - N_3 - N_2 + N_1")
    print("where N_d is the number of primitive characters with conductor N whose order divides d.")
    
    print("\nStep 3: Calculate N_d for d = 6, 3, 2, 1.")
    print("N_d is multiplicative, so N_d(N) = product over p^a || N of N_d(p^a).")

    results_d = {}
    divisors_of_6 = [6, 3, 2, 1]
    
    calc_strings = {}
    
    for d in divisors_of_6:
        total_Nd = 1
        current_calc_list = []
        for p, a in factors:
            count = count_prim_chars_order_divides_d(p, a, d)
            total_Nd *= count
            current_calc_list.append(f"N_{d}({p}^{a})={count}")
        results_d[d] = total_Nd
        calc_strings[d] = " * ".join(current_calc_list)

    n6, n3, n2, n1 = results_d[6], results_d[3], results_d[2], results_d[1]

    print("\nCalculating components N_d(p^a):")
    for d in divisors_of_6:
        print(f"For d={d}:")
        print(f"  N_{d}({N}) = {calc_strings[d]} = {results_d[d]}")
    
    final_result = n6 - n3 - n2 + n1

    print("\nStep 4: Apply the inclusion-exclusion formula.")
    print("Number of order 6 characters = N_6 - N_3 - N_2 + N_1")
    print(f"= {n6} - {n3} - {n2} + {n1}")
    print(f"= {final_result}")

solve_problem()
<<<100>>>