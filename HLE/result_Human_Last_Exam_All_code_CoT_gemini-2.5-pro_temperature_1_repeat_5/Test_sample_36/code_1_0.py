import math
from itertools import combinations, product

def get_divisors(n):
    """Returns the sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def lcm(a, b):
    """Returns the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def get_element_orders(cyclic_factors):
    """
    Calculates the set of possible element orders for a direct product of cyclic groups.
    E.g., for Z_6 x Z_3, cyclic_factors is [6, 3].
    """
    if not cyclic_factors:
        return {1}
    
    orders_per_factor = [get_divisors(n) for n in cyclic_factors]
    all_orders = set()
    
    # Generate all combinations of orders, one from each factor
    order_combinations = product(*orders_per_factor)
    
    for combo in order_combinations:
        current_lcm = 1
        for order in combo:
            current_lcm = lcm(current_lcm, order)
        all_orders.add(current_lcm)
        
    return sorted(list(all_orders))

def find_antichains(s):
    """Finds all antichains in a set s with respect to divisibility."""
    antichains_set = set()
    # Iterate through the power set of s (non-empty subsets)
    for i in range(1, len(s) + 1):
        for subset in combinations(s, i):
            is_antichain = True
            # Check if any element in the subset divides another
            if len(subset) > 1:
                for j, k in combinations(subset, 2):
                    if j % k == 0 or k % j == 0:
                        is_antichain = False
                        break
            if is_antichain:
                antichains_set.add(frozenset(subset))
    return antichains_set

# Step 1: Define the two non-isomorphic Abelian groups of order 18.
# G1 = Z_18
# G2 = Z_6 x Z_3

# Step 2: For each group, find the set of possible character orders.
orders_g1 = get_element_orders([18])
orders_g2 = get_element_orders([6, 3])

# Step 3: Find the set of all antichains for each set of orders.
antichains1 = find_antichains(orders_g1)
antichains2 = find_antichains(orders_g2)

# Step 4: Combine the results. The total number of unique sets is the size of the union
# of the sets of antichains.
all_unique_antichains = antichains1.union(antichains2)
total_unique_sets = len(all_unique_antichains)
num_antichains1 = len(antichains1)
num_antichains2 = len(antichains2)
num_common_antichains = len(antichains1.intersection(antichains2))

print(f"The set of character orders for G1 (Z_18) is: {orders_g1}")
print(f"The number of unique sets from G1 is the number of antichains: {num_antichains1}")
print(f"\nThe set of character orders for G2 (Z_6 x Z_3) is: {orders_g2}")
print(f"The number of unique sets from G2 is the number of antichains: {num_antichains2}")
print(f"\nNumber of sets also found for G1: {num_common_antichains}")

print("\nThe total number of unique sets is the size of the union of the sets from G1 and G2.")
print("Using the Principle of Inclusion-Exclusion, the total number is:")
print(f"{num_antichains1} + {num_antichains2} - {num_common_antichains} = {total_unique_sets}")
<<<9>>>