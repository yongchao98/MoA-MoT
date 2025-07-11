import math
from itertools import combinations

def get_character_orders_for_groups():
    """
    Determines the set of possible character orders for each non-isomorphic
    Abelian group of order 18. The dual group G_hat is isomorphic to G,
    so character orders are the same as element orders.
    """
    # Group 1: G1 = Z_18 (cyclic group)
    # The orders of elements (and characters) are the divisors of 18.
    n1 = 18
    orders_g1 = set()
    for i in range(1, int(math.sqrt(n1)) + 1):
        if n1 % i == 0:
            orders_g1.add(i)
            orders_g1.add(n1 // i)

    # Group 2: G2 = Z_2 x Z_3 x Z_3
    # The order of an element (a, b, c) is lcm(ord(a), ord(b), ord(c)).
    # ord(a) in {1, 2}, ord(b) in {1, 3}, ord(c) in {1, 3}.
    orders_g2 = set()
    for o1 in {1, 2}:
        for o2 in {1, 3}:
            for o3 in {1, 3}:
                lcm_val = (o1 * o2) // math.gcd(o1, o2)
                lcm_val = (lcm_val * o3) // math.gcd(lcm_val, o3)
                orders_g2.add(lcm_val)
                
    return sorted(list(orders_g1)), sorted(list(orders_g2))

def count_antichains(elements):
    """
    An antichain is a subset of elements where no element divides another.
    This function finds all such non-empty antichains.
    """
    antichains = set()
    num_elements = len(elements)
    # Iterate through all non-empty subsets of the elements
    for i in range(1, 1 << num_elements):
        subset = []
        for j in range(num_elements):
            if (i >> j) & 1:
                subset.append(elements[j])
        
        # Check if the current subset is an antichain
        is_ac = True
        if len(subset) > 1:
            for x, y in combinations(subset, 2):
                if x % y == 0 or y % x == 0:
                    is_ac = False
                    break
        
        if is_ac:
            antichains.add(frozenset(subset))
            
    return antichains

# --- Main execution ---

print("Step 1: Determine the character orders for Abelian groups of order 18.")
orders_g1, orders_g2 = get_character_orders_for_groups()
print(f"For G1 = Z_18, the possible character orders are: {orders_g1}")
print(f"For G2 = Z_2 x Z_3 x Z_3, the possible character orders are: {orders_g2}")
print("-" * 50)

print("Step 2: Count the number of unique eigenvalue sets for each group.")
print("This is equivalent to counting the number of antichains in the set of character orders (under divisibility).")

antichains_g1 = count_antichains(orders_g1)
num_sets_g1 = len(antichains_g1)
print(f"For G1, there are {num_sets_g1} unique sets.")

antichains_g2 = count_antichains(orders_g2)
num_sets_g2 = len(antichains_g2)
print(f"For G2, there are {num_sets_g2} unique sets.")
print("-" * 50)

print("Step 3: Calculate the total number of unique sets across all groups.")
# The total number is the size of the union of the sets of antichains
combined_antichains = antichains_g1.union(antichains_g2)
num_intersection = len(antichains_g1.intersection(antichains_g2))
total_unique_sets = len(combined_antichains)

print("The final count is the size of the union of the sets from G1 and G2.")
print("Using the Principle of Inclusion-Exclusion:")
print(f"|G1_sets U G2_sets| = |G1_sets| + |G2_sets| - |G1_sets intersect G2_sets|")
print(f"Total unique sets = {num_sets_g1} + {num_sets_g2} - {num_intersection}")
print(f"Total unique sets = {total_unique_sets}")

<<<9>>>