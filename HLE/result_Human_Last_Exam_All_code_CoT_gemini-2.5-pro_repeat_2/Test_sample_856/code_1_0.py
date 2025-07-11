import math

def get_prime_factorization(n):
    """Computes the prime factorization of n."""
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

def main():
    """
    Determines the number of closed orientable 3-manifolds with a fundamental group
    of cardinality 10!.
    """
    
    # Step 1 & 2: State the topological and group-theoretic conditions.
    print("Step 1: A closed orientable 3-manifold with a finite fundamental group G must be a spherical manifold.")
    print("Step 2: The group G must have the property that all of its abelian subgroups are cyclic.")
    
    # Step 3: State the consequence for Sylow p-subgroups.
    print("Step 3: This implies that for any odd prime p, the Sylow p-subgroups of G must be cyclic.")
    
    # Step 4: Calculate the order of the group and its prime factorization.
    group_order_val = math.factorial(10)
    print(f"\nStep 4: The order of the fundamental group is 10! = {group_order_val}.")
    prime_factors = get_prime_factorization(group_order_val)
    factors_str = " * ".join([f"{p}^{e}" for p, e in prime_factors.items()])
    print(f"The prime factorization of {group_order_val} is {factors_str}.")
    
    # Step 5: Focus on the Sylow 3-subgroups.
    p = 3
    sylow_3_order = prime_factors[p]**p # whoops, this is wrong. It's p^exponent
    sylow_3_order = p**prime_factors[p]
    print(f"\nStep 5: For the prime p = {p}, the order of a Sylow p-subgroup is {p}^{prime_factors[p]} = {sylow_3_order}.")
    print(f"For G to be a valid manifold group, its Sylow {p}-subgroups must be cyclic.")

    # Step 6 & 7: Explain the contradiction.
    print("\nStep 6: A cyclic group of order 81 must contain an element of order 81.")
    print("Step 7: However, it is a fact of finite group theory that NO group of order 10! can have an element of order 81.")
    print("Therefore, for any group G of order 10!, its Sylow 3-subgroups cannot be cyclic.")

    # Step 8 & 9: Draw the final conclusion.
    print("\nStep 8: A non-cyclic p-group (for p odd) is known to contain a subgroup isomorphic to Z_p x Z_p.")
    print(f"In our case, any group of order 10! must contain a subgroup isomorphic to Z_{p} x Z_{p} = Z_3 x Z_3.")
    print("Step 9: The subgroup Z_3 x Z_3 is abelian but not cyclic, which violates the condition from Step 2.")
    
    # Final Answer
    num_manifolds = 0
    print("\nConclusion: No group of order 10! can be the fundamental group of a closed orientable 3-manifold.")
    print(f"The number of such manifolds is {num_manifolds}.")
    
if __name__ == "__main__":
    main()
