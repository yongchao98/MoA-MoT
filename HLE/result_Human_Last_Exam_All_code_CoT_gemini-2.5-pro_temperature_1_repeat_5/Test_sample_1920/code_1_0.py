import math

def get_prime_factorization(n):
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
    """Calculates Euler's totient function phi(n)."""
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

def solve():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Let d = {d}. We want to find the number of primitive Dirichlet characters of conductor d and order {g}.")
    print("-" * 80)

    # Step 1: Prime factorization of d
    print(f"Step 1: Find the prime factorization of d = {d}.")
    d_factors = get_prime_factorization(d)
    primes = sorted(list(d_factors.keys()))
    print(f"The prime factorization of {d} is: {' * '.join(map(str, primes))}.")
    print("\nSince d is square-free, a character chi modulo d is a product of characters chi_p for each prime p dividing d.")
    print("For chi to be primitive with conductor d, each chi_p must be non-principal (order > 1).")
    print(f"The order of chi is lcm(ord(chi_p) for p|d), which must be {g}.")
    print("-" * 80)

    # Step 2: Check character group orders
    print("Step 2: Check the orders of the character groups modulo each prime factor.")
    phi_values = {p: p - 1 for p in primes}
    print("The group of characters modulo a prime p is cyclic of order phi(p) = p-1.")
    for p in primes:
        print(f"  - For p = {p}, the character group has order phi({p}) = {phi_values[p]}.")
    print("-" * 80)

    # Step 3: Determine possible orders for component characters
    print("Step 3: Determine the possible orders for each component character chi_p.")
    print(f"The order of each chi_p must be greater than 1 and be a divisor of {g}=6.")
    possible_orders = [k for k in [2, 3, 6] if g % k == 0]
    print(f"So, the possible orders for each chi_p are {possible_orders}.")
    
    # Check if these orders are valid
    for p in primes:
        for k in possible_orders:
            if phi_values[p] % k != 0:
                # This case is not hit for the given problem but is good practice.
                print(f"Error: For p={p}, order {k} is not possible as it does not divide phi({p})={phi_values[p]}.")
                return
    print(f"For each prime factor p, phi(p) is divisible by all k in {possible_orders}.")
    print("-" * 80)
    
    # Step 4: Count characters for each possible order
    print("Step 4: Count the number of characters for each allowed order.")
    print("The number of characters of order k modulo p is phi(k).")
    phi_k_values = {k: phi(k) for k in possible_orders}
    for k, val in phi_k_values.items():
        print(f"  - Number of characters of order {k}: phi({k}) = {val}")

    total_per_prime = sum(phi_k_values.values())
    print(f"\nFor each prime factor, there are { ' + '.join(map(str, phi_k_values.values())) } = {total_per_prime} non-principal characters whose order divides {g}.")
    print("-" * 80)
    
    # Step 5: Apply inclusion-exclusion principle
    print("Step 5: Use the Principle of Inclusion-Exclusion to find the number of characters with order exactly 6.")
    num_primes = len(primes)
    
    # Total primitive characters with order dividing 6
    total_combinations = total_per_prime ** num_primes
    print(f"The total number of primitive characters with order dividing 6 is the product of the counts for each prime:")
    print(f"This is {total_per_prime}^{num_primes} = {total_combinations}.")
    
    print("\nFrom this, we subtract the cases where the lcm of orders is a proper divisor of 6 (i.e., 2 or 3).")
    
    # Cases where lcm divides 2 (i.e., all orders are 2)
    num_ord_2 = phi_k_values[2]
    sub_lcm_2 = num_ord_2 ** num_primes
    print(f"  - Number of combinations where order is 2: (phi(2))^{num_primes} = {num_ord_2}^{num_primes} = {sub_lcm_2}.")
    
    # Cases where lcm divides 3 (i.e., all orders are 3)
    num_ord_3 = phi_k_values[3]
    sub_lcm_3 = num_ord_3 ** num_primes
    print(f"  - Number of combinations where order is 3: (phi(3))^{num_primes} = {num_ord_3}^{num_primes} = {sub_lcm_3}.")

    final_count = total_combinations - sub_lcm_2 - sub_lcm_3
    print("\nThe final count is obtained by subtracting these from the total.")
    print(f"Final Answer = (Total combinations) - (Combinations with order 2) - (Combinations with order 3)")
    print(f"             = {total_combinations} - {sub_lcm_2} - {sub_lcm_3}")
    print(f"             = {final_count}")
    print("-" * 80)
    
    return final_count

# Run the solver
final_answer = solve()
print(f"<<<{final_answer}>>>")
