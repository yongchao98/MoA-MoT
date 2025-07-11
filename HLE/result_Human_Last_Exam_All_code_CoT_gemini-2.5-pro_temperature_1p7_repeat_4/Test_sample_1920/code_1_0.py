import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary.
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
    Computes Euler's totient function phi(n).
    """
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

def solve():
    """
    Calculates the number of primitive Dirichlet characters of conductor d and order g.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.\n")

    # Step 1: Factor the conductor
    factors = get_prime_factorization(d)
    primes = list(factors.keys())
    factor_str = " * ".join(map(str, primes))
    print(f"Step 1: The prime factorization of d = {d} is {factor_str}.")

    # Check if square-free
    is_square_free = all(v == 1 for v in factors.values())
    if not is_square_free:
        print("Conductor d is not square-free. The calculation is more complex. This script handles square-free conductors.")
        return

    print("The conductor d is square-free. This simplifies the problem.")
    print("A character chi is primitive iff its components chi_p are all non-principal (order > 1).\n")

    # Step 2: Analyze order conditions
    print("Step 2: Analyzing the order conditions for component characters.")
    print(f"The order of chi is lcm(ord(chi_p) for p in {primes}). We need this lcm to be {g}.")
    
    order_divisors = [k for k in range(1, g + 1) if g % k == 0]
    primitive_orders = [k for k in order_divisors if k > 1]
    
    print(f"This implies ord(chi_p) must divide {g}. For primitivity, ord(chi_p) must be > 1.")
    print(f"So, the possible orders for each component character are {primitive_orders}.\n")

    # Step 3: Check validity and count characters for each order
    print("Step 3: Counting the number of choices for each component character.")
    valid = True
    for p in primes:
        phi_p = p - 1
        for k in primitive_orders:
            if phi_p % k != 0:
                print(f"Error: Order {k} is not a divisor of phi({p}) = {phi_p}.")
                valid = False
    if not valid:
        return
        
    print(f"For each prime factor p, the orders {primitive_orders} are valid since they divide phi(p)=p-1.")

    phi_values = {k: phi(k) for k in primitive_orders}
    print("The number of characters of order k is phi(k):")
    for k, v in phi_values.items():
        print(f"  phi({k}) = {v}")

    # Step 4: Use inclusion-exclusion to find the final count
    print("\nStep 4: Using combinatorial counting to find the total number.")

    # Total primitive characters with order dividing g=6
    # This means ord(chi_p) is in {2, 3, 6} for each component
    num_choices_per_component_div_6 = sum(phi_values.values())
    total_dividing_6 = num_choices_per_component_div_6 ** len(primes)
    
    # Number of primitive characters with order dividing 2 (must be 2)
    # This means ord(chi_p) = 2 for each component
    num_choices_order_2 = phi_values.get(2, 0)
    total_order_2 = num_choices_order_2 ** len(primes)

    # Number of primitive characters with order dividing 3 (must be 3)
    # This means ord(chi_p) = 3 for each component
    num_choices_order_3 = phi_values.get(3, 0)
    total_order_3 = num_choices_order_3 ** len(primes)

    # The result is by principle of inclusion-exclusion
    result = total_dividing_6 - total_order_2 - total_order_3
    
    print("Let N_k be the number of primitive characters of order k.")
    print(f"Number of primitive characters with order dividing 6 (N_2 + N_3 + N_6) = ({' + '.join(map(str, phi_values.values()))})^{len(primes)} = {num_choices_per_component_div_6}^{len(primes)} = {total_dividing_6}.")
    print(f"Number of primitive characters of order 2 (N_2) = {num_choices_order_2}^{len(primes)} = {total_order_2}.")
    print(f"Number of primitive characters of order 3 (N_3) = {num_choices_order_3}^{len(primes)} = {total_order_3}.")
    
    print("\nThe number of characters of order 6 is (Total with order dividing 6) - N_2 - N_3.")
    print(f"Final Calculation: {total_dividing_6} - {total_order_2} - {total_order_3} = {result}")

    print("\nThe final equation is:")
    phi2, phi3, phi6 = phi_values[2], phi_values[3], phi_values[6]
    n_primes = len(primes)
    print(f"(phi(2) + phi(3) + phi(6))^{n_primes} - (phi(2))^{n_primes} - (phi(3))^{n_primes}")
    print(f"= ({phi2} + {phi3} + {phi6})^{n_primes} - {phi2}^{n_primes} - {phi3}^{n_primes}")
    print(f"= {num_choices_per_component_div_6}^{n_primes} - {total_order_2} - {total_order_3}")
    print(f"= {total_dividing_6} - {total_order_2} - {total_order_3} = {result}\n")
    print(f"The number of primitive Dirichlet characters of conductor {d} and order {g} is {result}.")

    return result

# Execute the solver
final_answer = solve()
print(f'<<<{final_answer}>>>')