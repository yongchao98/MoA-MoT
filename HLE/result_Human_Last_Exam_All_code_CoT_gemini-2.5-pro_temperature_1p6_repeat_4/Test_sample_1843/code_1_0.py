import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {p: a}.
    """
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

def euler_phi_from_factors(factors):
    """
    Computes Euler's totient function from prime factors.
    """
    result = 1
    for p, a in factors.items():
        result *= (p - 1) * (p**(a - 1))
    return result
    
def euler_phi(n):
    """
    Computes Euler's totient function.
    """
    if n == 1:
        return 1
    return euler_phi_from_factors(get_prime_factorization(n))

def get_divisors(n):
    """
    Returns a sorted list of divisors of n.
    """
    divs = {1, n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve():
    """
    Solves the problem of counting primitive Dirichlet characters.
    """
    N = 36036
    target_order = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order {target_order}.")
    print("-" * 30)

    # 1. Prime factorization of N
    factors_N = get_prime_factorization(N)
    conductors = [p**a for p, a in factors_N.items()]
    print(f"The conductor N = {N} factorizes into prime powers: { {p: a for p, a in factors_N.items()} }.")
    print(f"The primitive component characters will have conductors: {conductors}\n")

    # 2. Count characters for each conductor
    print("Step 1: For each conductor, count the number of primitive characters whose order divides 6.\n")

    total_chars = 1
    char_counts = {}

    # Conductor d = 4 (2^2)
    d = 4
    phi_d = euler_phi(d)
    primitive_orders = [k for k in get_divisors(phi_d) if k > 1] # Primitive orders mod 4 are {2}
    allowed_orders = [k for k in primitive_orders if target_order % k == 0]
    count = sum(euler_phi(k) for k in allowed_orders)
    char_counts[d] = count
    total_chars *= count
    print(f"For conductor d = {d}:")
    print(f"  The group of characters has order phi({d}) = {phi_d}.")
    print(f"  The primitive characters have orders: {primitive_orders}.")
    print(f"  Orders that divide {target_order}: {allowed_orders}.")
    print(f"  Number of characters = phi(2) = {count}.\n")

    # Conductor d = 9 (3^2)
    d = 9
    phi_d = euler_phi(d)
    phi_d_prev = euler_phi(3)
    primitive_orders = [k for k in get_divisors(phi_d) if phi_d_prev % k != 0] # Orders divide 6 but not 2
    allowed_orders = [k for k in primitive_orders if target_order % k == 0]
    count = sum(euler_phi(k) for k in allowed_orders)
    char_counts[d] = count
    total_chars *= count
    print(f"For conductor d = {d}:")
    print(f"  The group of characters has order phi({d}) = {phi_d}.")
    print(f"  The primitive characters have orders dividing phi(9) but not phi(3): {primitive_orders}.")
    print(f"  Orders that divide {target_order}: {allowed_orders}.")
    print(f"  Number of characters = phi(3) + phi(6) = {euler_phi(3)} + {euler_phi(6)} = {count}.\n")

    # Conductors d = p (prime)
    for p in [7, 11, 13]:
        d = p
        phi_d = euler_phi(d)
        primitive_orders = [k for k in get_divisors(phi_d) if k > 1]
        allowed_orders = [k for k in primitive_orders if target_order % k == 0]
        count = sum(euler_phi(k) for k in allowed_orders)
        char_counts[d] = count
        total_chars *= count
        print(f"For conductor d = {d}:")
        print(f"  The group of characters has order phi({d}) = {phi_d}.")
        print(f"  The primitive (non-trivial) characters have orders: {primitive_orders}.")
        print(f"  Orders that divide {target_order}: {allowed_orders}.")
        
        counts_str = " + ".join([f"phi({k})" for k in allowed_orders])
        values_str = " + ".join([str(euler_phi(k)) for k in allowed_orders])
        if len(allowed_orders) > 1:
            print(f"  Number of characters = {counts_str} = {values_str} = {count}.\n")
        else:
            print(f"  Number of characters = {counts_str} = {count}.\n")

    # 3. Final calculation logic
    print("Step 2: Check if the resulting character order must be 6.\n")
    print("Let chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13 be a character formed by primitive components.")
    print("The order of chi is lcm(ord(chi_4), ord(chi_9), ...).")
    print(" - For conductor 4, the order of chi_4 must be 2. So, ord(chi) is divisible by 2.")
    print(" - For conductor 9, the order of chi_9 must be 3 or 6. So, ord(chi) is divisible by 3.")
    print("Therefore, the order of any such character must be divisible by lcm(2, 3) = 6.")
    print("Since all component orders divide 6, their lcm must also divide 6.")
    print("A number that is a multiple of 6 and a divisor of 6 must be exactly 6.")
    print("So, we just need to multiply the counts for each conductor.\n")
    
    print("Step 3: Final Calculation\n")
    # 4. Final result
    calc_str = " * ".join(map(str, char_counts.values()))
    print(f"The total number of primitive characters of order 6 is the product of the individual counts:")
    print(f"Total = {char_counts[4]} * {char_counts[9]} * {char_counts[7]} * {char_counts[11]} * {char_counts[13]} = {total_chars}")

solve()

<<<100>>>