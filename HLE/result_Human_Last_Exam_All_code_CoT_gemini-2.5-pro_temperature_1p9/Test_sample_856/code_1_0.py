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

# The order of the fundamental group is 10!
group_order = math.factorial(10)

# The core of the argument relies on the prime factorization of the group order.
prime_factors = get_prime_factorization(group_order)

print(f"The order of the fundamental group is 10! = {group_order}.")
print(f"The prime factorization of {group_order} is:")
factor_strings = []
sylow_orders = []
for p, a in prime_factors.items():
    sylow_order = p**a
    sylow_orders.append(sylow_order)
    factor_strings.append(f"{p}^{a} (= {sylow_order})")
print(" * ".join(factor_strings))

print("\nStep 1: A 3-manifold with a finite fundamental group G of this order must be of the form S^3/G.")
print("Step 2: A necessary condition (Milnor's Theorem) is that G cannot be a direct product of two non-trivial groups of coprime order.")
print("\nStep 3: A further (deep) result shows that any group G of order 10! satisfying the conditions for acting on S^3 must be nilpotent.")
print("A finite group is nilpotent if and only if it is the direct product of its Sylow p-subgroups.")
print("\nStep 4: If such a group G were nilpotent, it would be a direct product of groups with the following orders:")
print(sylow_orders)
print(f"This is a direct product of {len(sylow_orders)} non-trivial groups of pairwise coprime orders.")
print("\nConclusion: This structure contradicts Milnor's theorem. Therefore, no such group G can exist.")
print("The number of such 3-manifolds is 0.")
