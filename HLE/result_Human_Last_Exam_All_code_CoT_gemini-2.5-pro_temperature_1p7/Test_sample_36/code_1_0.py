import math
from itertools import chain, combinations

print("This script calculates the number of unique eigenvalue sets on the unit circle for representations of Abelian groups of order 18.")

# Step 1: Explain the theoretical background.
print("\nStep 1: For any representation of a finite group, all eigenvalues are roots of unity and lie on the unit circle D.")
print("Thus, S(rho) intersect D is simply S(rho), the set of all eigenvalues.")

print("\nStep 2: For an Abelian group, S(rho) is a union of the images of its constituent characters.")
print("The image of a character is a set of k-th roots of unity, U_k, where k is the character's order.")

# Step 2: Determine possible character orders for each group type.
print("\nStep 3: The Abelian groups of order 18 are G1 = Z_18 and G2 = Z_6 x Z_3.")
# For G1 = Z_18, orders are divisors of 18.
orders_g1 = {d for d in range(1, 19) if 18 % d == 0}
print(f"The possible character orders for G1 = Z_18 are the divisors of 18: {sorted(list(orders_g1))}")

# For G2 = Z_6 x Z_3, orders are lcm of element orders from Z_6 and Z_3.
# The element orders in Z_6 are {1,2,3,6} and in Z_3 are {1,3}.
orders_g2 = set()
for o1 in {1, 2, 3, 6}:
    for o2 in {1, 3}:
        orders_g2.add((o1 * o2) // math.gcd(o1, o2))
print(f"The possible character orders for G2 = Z_6 x Z_3 are: {sorted(list(orders_g2))}")

# Step 3: Combine orders to get the full set for the poset.
poset_elements = sorted(list(orders_g1.union(orders_g2)))
print(f"\nStep 4: The total pool of possible orders k for image sets U_k is the union: P = {poset_elements}")
print("The unique eigenvalue sets correspond to the non-empty antichains of the poset (P, |), where '|' is divisibility.")

# Step 4: Find and count all non-empty antichains.
all_subsets = chain.from_iterable(combinations(poset_elements, r) for r in range(1, len(poset_elements) + 1))

antichains_by_size = {}

for subset in all_subsets:
    is_antichain = True
    if len(subset) > 1:
        # Check all pairs for comparability (divisibility).
        for e1, e2 in combinations(subset, 2):
            if (e1 % e2 == 0) or (e2 % e1 == 0):
                is_antichain = False
                break
    
    if is_antichain:
        size = len(subset)
        if size not in antichains_by_size:
            antichains_by_size[size] = []
        antichains_by_size[size].append(subset)

print("\nStep 5: Counting the non-empty antichains by their size.")
total_count = 0
equation_parts = []
for size in sorted(antichains_by_size.keys()):
    count = len(antichains_by_size[size])
    print(f"Number of antichains of size {size}: {count}")
    total_count += count
    equation_parts.append(str(count))

print("\nFinal Result:")
print("The total number of unique sets is the sum of the counts of antichains of each size.")
final_equation = " + ".join(equation_parts)
print(f"Final Equation: {final_equation} = {total_count}")

# Final Answer Block
print("\n<<<9>>>")