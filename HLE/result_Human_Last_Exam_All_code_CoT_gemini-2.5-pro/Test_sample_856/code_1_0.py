import math

def get_prime_power(n, p):
    """Calculates the exponent of prime p in the factorization of n! using Legendre's Formula."""
    power = 0
    i = 1
    while True:
        term = math.floor(n / (p**i))
        if term == 0:
            break
        power += term
        i += 1
    return power

def main():
    """
    Solves the problem by analyzing the properties of the fundamental group.
    """
    n = 10
    order_val = math.factorial(n)
    
    # Step 1: Explain the topological condition
    print("Step 1: The Topological Condition")
    print("A closed orientable 3-manifold with a finite fundamental group, Gamma, must be a spherical space form.")
    print("This means Gamma must be a finite group that can act freely on the 3-sphere S^3.")
    print("-" * 20)
    
    # Step 2: State the group-theoretic condition
    print("Step 2: The Group-Theoretic Condition")
    print("A group Gamma can act freely on S^3 only if it has no subgroups isomorphic to C_p x C_p for any prime p.")
    print("(C_p is the cyclic group of order p).")
    print("-" * 20)

    # Step 3: Analyze the group order
    print("Step 3: Analyze the Group of Order 10!")
    print(f"The order of the fundamental group Gamma is given as {n}! = {order_val}.")
    
    # Step 4: Focus on prime p=3
    p = 3
    power_of_3 = get_prime_power(n, p)
    sylow_3_order = p**power_of_3
    
    print(f"Let's check the condition for the prime p = {p}.")
    print(f"The order of a Sylow {p}-subgroup of a group of order {n}! is {p}^{power_of_3} = {sylow_3_order}.")
    print("A known theorem states that if a p-group (for odd p) is not cyclic, it must contain a subgroup isomorphic to C_p x C_p.")
    print(f"Therefore, for Gamma to be a valid fundamental group, its Sylow {p}-subgroups (of order {sylow_3_order}) must be cyclic.")
    print("-" * 20)
    
    # Step 5: Check the condition for S_10
    print(f"Step 4: Check the Condition for the Symmetric Group S_{n}")
    print(f"A group with a cyclic subgroup of order {sylow_3_order} must contain an element of order {sylow_3_order}.")
    print(f"Let's check if the symmetric group S_{n} (which has order {n}!) contains an element of order {sylow_3_order}.")
    print(f"The order of an element in S_{n} is the least common multiple (lcm) of its disjoint cycle lengths, which must sum to {n}.")
    print(f"To have an element of order {sylow_3_order}, we would need a cycle whose length is a multiple of {sylow_3_order}. This is impossible in S_{n} since the maximum cycle length is {n}.")
    print(f"Therefore, S_{n} has no element of order {sylow_3_order}.")
    print("-" * 20)
    
    # Step 6: Conclusion
    print("Step 5: Conclusion")
    print(f"Since S_{n} has no element of order {sylow_3_order}, its Sylow {p}-subgroups are not cyclic. This means S_{n} contains a C_{p} x C_{p} subgroup and cannot be the fundamental group.")
    print("This property (of having non-cyclic Sylow subgroups for relevant primes) holds for any group of order 10!.")
    print("Therefore, no group of order 10! can be the fundamental group of a closed orientable 3-manifold.")
    print("The number of such manifolds is 0.")

if __name__ == "__main__":
    main()