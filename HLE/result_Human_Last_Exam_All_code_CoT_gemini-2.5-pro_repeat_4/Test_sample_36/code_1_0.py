import math
from itertools import chain, combinations

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def divides(a, b):
    """Checks if a divides b, where a is not equal to b."""
    return a != b and b % a == 0

def get_antichains(S):
    """
    Finds all antichains in a set S with respect to the "divides" relation.
    An antichain is a subset where no element divides another.
    """
    antichains = set()
    # Generate all non-empty subsets of S
    all_subsets = chain.from_iterable(combinations(S, r) for r in range(1, len(S) + 1))
    
    for subset in all_subsets:
        is_ac = True
        s_list = list(subset)
        # Check if any element in the subset divides another
        for i in range(len(s_list)):
            for j in range(i + 1, len(s_list)):
                a, b = s_list[i], s_list[j]
                if divides(a, b) or divides(b, a):
                    is_ac = False
                    break
            if not is_ac:
                break
        
        if is_ac:
            # Add the antichain as a frozenset to make it hashable for the set
            antichains.add(frozenset(subset))
            
    return antichains

# Step 1: Identify the character orders for the two Abelian groups of order 18.
# Group 1: G1 = Z_18. Character orders are the divisors of 18.
char_orders_g1 = get_divisors(18)

# Group 2: G2 = Z_2 x Z_3 x Z_3. Exponent is lcm(2, 3, 3) = 6.
# Character orders are the divisors of the exponent: {1, 2, 3, 6}.
char_orders_g2 = get_divisors(6)

# Step 2: Find the number of unique sets S(rho) for each group by counting antichains.
antichains_g1 = get_antichains(char_orders_g1)
num_sets_g1 = len(antichains_g1)

antichains_g2 = get_antichains(char_orders_g2)
num_sets_g2 = len(antichains_g2)

# Step 3: Calculate the total number of unique sets across both groups.
# This is the size of the union of the two sets of antichains.
total_unique_sets = antichains_g1.union(antichains_g2)
num_total_unique_sets = len(total_unique_sets)

# The number of common sets is the size of the intersection.
num_common_sets = len(antichains_g1.intersection(antichains_g2))

# Step 4: Print the final calculation.
print(f"Number of unique sets S(rho) for Z_18: {num_sets_g1}")
print(f"Number of unique sets S(rho) for Z_2 x Z_3 x Z_3: {num_sets_g2}")
print(f"Number of common sets between the two groups: {num_common_sets}")
print("The total number of unique sets is given by the Principle of Inclusion-Exclusion:")
print(f"Total = (Sets from Z_18) + (Sets from Z_2 x Z_3 x Z_3) - (Common Sets)")
print(f"Total = {num_sets_g1} + {num_sets_g2} - {num_common_sets} = {num_total_unique_sets}")
<<<9>>>