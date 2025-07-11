import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary {prime: exponent}."""
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

# Step 1: Define the order of the fundamental group.
n = 10
order = math.factorial(n)
print(f"The order of the fundamental group is {n}! = {order}")

# Step 2: Find the prime factorization of the order.
prime_factors = get_prime_factorization(order)
factorization_str = " * ".join([f"{p}^{e}" for p, e in sorted(prime_factors.items())])
print(f"The prime factorization of {order} is: {factorization_str}")
print("-" * 30)

# Step 3: Analyze the Sylow subgroups based on the group-theoretic conditions
# for a group to be the fundamental group of a spherical 3-manifold.
# - Sylow p-subgroups for odd p must be cyclic.
# - Sylow 2-subgroup must be cyclic or generalized quaternion.

possibilities_per_prime = {}
print("Analyzing the required structure of the Sylow subgroups:")

for p, e in sorted(prime_factors.items()):
    sylow_order = p**e
    print(f"\nFor prime p = {p}:")
    print(f"  The Sylow {p}-subgroup has order {p}^{e} = {sylow_order}.")
    
    if p % 2 != 0:  # Odd prime
        # For an odd prime, the Sylow p-subgroup must be cyclic.
        # There is exactly one cyclic group of order p^e up to isomorphism.
        possibilities_per_prime[p] = 1
        print(f"  The condition requires this subgroup to be cyclic. This gives 1 possibility.")
    else:  # p = 2
        # For p=2, the Sylow 2-subgroup can be cyclic or generalized quaternion.
        # A generalized quaternion group of order 2^e exists if e >= 3.
        if e >= 3:
            possibilities_per_prime[p] = 2
            print(f"  The condition requires this subgroup to be cyclic (C_{sylow_order}) or generalized quaternion (Q_{sylow_order}).")
            print(f"  This gives 2 possibilities.")
        else:
            # If e < 3, a generalized quaternion group of order 2^e does not exist.
            # The only option is cyclic.
            possibilities_per_prime[p] = 1
            print(f"  The exponent {e} < 3, so the subgroup must be cyclic. This gives 1 possibility.")

print("-" * 30)

# Step 4: Calculate the total number of non-isomorphic groups.
# A group with these properties is a direct product of its Sylow subgroups.
# The total number of groups is the product of possibilities for each Sylow subgroup.
total_manifolds = 1
for p in sorted(possibilities_per_prime.keys()):
    total_manifolds *= possibilities_per_prime[p]
    
equation_parts = [str(v) for k, v in sorted(possibilities_per_prime.items())]
equation_str = " * ".join(equation_parts)

print("The total number of manifolds is the product of the possibilities for each prime's Sylow subgroup.")
print(f"Final calculation: {equation_str} = {total_manifolds}")
print(f"\nTherefore, there are {total_manifolds} closed orientable 3-manifolds with a fundamental group of order 10!.")
