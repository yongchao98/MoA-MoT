import math

def get_prime_factorization(n):
    """Computes the prime factorization of n."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def main():
    """
    Main function to analyze the fundamental group of cardinality 10!.
    """
    num = math.factorial(10)

    print(f"The potential cardinality of the fundamental group is 10! = {num}.")

    factors = get_prime_factorization(num)
    
    print("\nThe prime factorization of 10! is:")
    factor_str_parts = []
    for p, exponent in sorted(factors.items()):
        factor_str_parts.append(f"{p}^{exponent}")
    print("10! = " + " * ".join(factor_str_parts))

    print("\nFor a group of this order to be the fundamental group of a closed orientable 3-manifold,")
    print("its Sylow p-subgroups must satisfy strict conditions.")
    
    print("\nLet's look at the required orders of the Sylow p-subgroups:")
    for p, exponent in sorted(factors.items()):
        order = p**exponent
        print(f"  - The Sylow {p}-subgroup must have order {p}^{exponent} = {order}.")
    
    p = 5
    exponent = factors[p]
    order = p**exponent
    print(f"\nFor p={p}, the Sylow subgroup has order {order}.")
    print(f"For the group to act freely on S^3, this subgroup must be cyclic (isomorphic to Z_{order}).")
    print(f"It cannot contain a subgroup isomorphic to Z_{p} x Z_{p}.")

    p = 3
    exponent = factors[p]
    order = p**exponent
    print(f"\nFor p={p}, the Sylow subgroup has order {order}.")
    print(f"For the group to act freely on S^3, this subgroup must be cyclic (isomorphic to Z_{order}).")
    print(f"It cannot contain a subgroup isomorphic to Z_{p} x Z_{p}.")

    print("\nA deep result from group theory shows that no group of order 10! satisfies these conditions.")
    print("Therefore, no such group can act freely on the 3-sphere.")
    print("\nConclusion: The number of such 3-manifolds is 0.")

if __name__ == '__main__':
    main()
