import math

def solve_manifold_problem():
    """
    Calculates the number of closed orientable 3-manifolds with a fundamental
    group of cardinality 10! by walking through the equivalent group theory problem.
    """
    
    # Step 1: Define the order of the fundamental group and find its prime factorization.
    n_factorial = 10
    group_order = math.factorial(n_factorial)

    def get_prime_factorization(num):
        factors = {}
        d = 2
        temp = num
        while d * d <= temp:
            while temp % d == 0:
                factors[d] = factors.get(d, 0) + 1
                temp //= d
            d += 1
        if temp > 1:
            factors[temp] = factors.get(temp, 0) + 1
        return factors

    prime_factors = get_prime_factorization(group_order)

    # Step 2: Print the reasoning based on the plan.
    print("This problem can be solved by translating it from topology to group theory.")
    print(f"The order of the fundamental group is 10! = {group_order}.")
    print("\n--- Step-by-step Reasoning ---")
    print("1. According to the Geometrization Theorem, a closed, orientable 3-manifold with a finite fundamental group G is unique up to homeomorphism. Thus, the problem is equivalent to counting the number of non-isomorphic groups of order 10! that can be such a fundamental group.")
    
    print("\n2. A finite group G can be a 3-manifold fundamental group if and only if it satisfies the 'pq-conditions': subgroups of order p^2 or pq (for primes p,q) must be cyclic.")

    print("\n3. A group G of order 10! satisfying these conditions must be a direct product of its Sylow subgroups. The prime factorization of 10! is:")
    print(f"   10! = {' * '.join([f'{p}^{e}' for p, e in prime_factors.items()])}")
    print("   So, G must be of the form: P_2 x P_3 x P_5 x P_7, where P_p is a Sylow p-subgroup.")

    print("\n4. We now count the number of valid choices for each Sylow subgroup:")
    
    choices = {}
    # For odd primes
    for p in sorted(prime_factors.keys()):
        if p % 2 != 0:
            order_p = p**prime_factors[p]
            print(f"   - For p={p}, the Sylow subgroup P_{p} has order {p}^{prime_factors[p]} = {order_p}.")
            print(f"     A p-group (for odd p) satisfying the conditions must be cyclic. There is only 1 choice (Z_{order_p}).")
            choices[p] = 1

    # For p=2
    p = 2
    order_2 = p**prime_factors[p]
    print(f"   - For p={p}, the Sylow subgroup P_{p} has order {p}^{prime_factors[p]} = {order_2}.")
    print(f"     A 2-group satisfying the conditions must be cyclic (Z_{order_2}) or generalized quaternion (Q_{order_2}).")
    print("     Thus, there are 2 choices for P_2.")
    choices[p] = 2
    
    # Step 5: Calculate the total number of groups and print the final equation.
    total_choices = 1
    for p in sorted(choices.keys()):
        total_choices *= choices[p]

    print("\n5. The total number of non-isomorphic groups is the product of the number of choices for each Sylow subgroup.")
    
    # The prompt requests that each number in the final equation be output.
    equation_parts = []
    for p in sorted(choices.keys()):
        equation_parts.append(f"{choices[p]}")

    print("Final Equation: " + " * ".join(equation_parts) + f" = {total_choices}")

    print(f"\nConclusion: There are {total_choices} such non-isomorphic groups, and therefore {total_choices} such 3-manifolds.")

solve_manifold_problem()