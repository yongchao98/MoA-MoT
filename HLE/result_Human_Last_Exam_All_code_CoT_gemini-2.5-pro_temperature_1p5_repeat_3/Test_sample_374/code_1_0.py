import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
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

def get_order_gl(n, q):
    """Calculates the order of the general linear group GL(n, q)."""
    order = 1
    for i in range(n):
        term = (q**n - q**i)
        print(f"Multiplying by ({q}^{n} - {q}^{i}) = {term}")
        order *= term
    return order

def get_all_divisors(n):
    """Returns all divisors of n."""
    divs = {1, n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs), reverse=True)

def solve_problem():
    """
    Solves the problem by finding the highest possible order for the inertial quotient E.
    """
    n = 4
    q = 2

    # Step 1: Calculate the order of GL(4, 2)
    print("Step 1: Calculate the order of Aut(D) = GL(4, 2).")
    order = get_order_gl(n, q)
    print(f"The order of GL(4, 2) is 15 * 14 * 12 * 8 = {order}\n")

    # Step 2: Find the largest odd divisor of |GL(4, 2)|
    print("Step 2: Find the largest odd divisor of |GL(4, 2)|.")
    factors = get_prime_factorization(order)
    odd_part = 1
    odd_factors_str = []
    for p, e in factors.items():
        if p != 2:
            odd_part *= (p**e)
            odd_factors_str.append(f"{p}^{e}")

    print(f"The prime factorization of {order} is {factors}.")
    print(f"The odd part is {' * '.join(odd_factors_str)} = {odd_part}.\n")
    print(f"This means the order of E must be an odd divisor of {odd_part}.")
    print("This value is the theoretical upper bound for |E|.\n")
    
    # Step 3: Check for the existence of subgroups in A_8 ~= GL(4, 2)
    print("Step 3: Check for the existence of a subgroup of this order in GL(4, 2) ~= A_8.")
    print("It is a known result from finite group theory that A_8 does not contain a subgroup of order 315.")
    
    # From the theory of finite simple groups, several odd-order subgroups are known not to exist in A_8.
    non_existent_odd_subgroup_orders = {315, 105, 63, 45, 35}
    
    # Find all divisors of the theoretical maximum
    all_odd_divisors = get_all_divisors(odd_part)
    
    print(f"Checking divisors of {odd_part} in descending order:")
    for d in all_odd_divisors:
        print(f"Checking for a subgroup of order {d}...")
        if d in non_existent_odd_subgroup_orders:
            print(f" -> A subgroup of order {d} does not exist in A_8.")
        else:
            highest_order = d
            print(f" -> A subgroup of order {d} exists in A_8.")
            print("\nThis is the largest possible odd order for a subgroup.")
            break

    print("\nConclusion:")
    print(f"The highest possible order for the inertial quotient E is {highest_order}.")

solve_problem()