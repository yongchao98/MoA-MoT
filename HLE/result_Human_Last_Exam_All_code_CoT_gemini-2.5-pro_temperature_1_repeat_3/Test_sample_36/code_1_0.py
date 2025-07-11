import numpy as np
from itertools import chain, combinations

def get_roots(d):
    """
    Returns a frozenset of the d-th roots of unity.
    We round the complex numbers to handle potential floating-point inaccuracies.
    """
    return frozenset(np.round(np.exp(2j * np.pi * k / d), 12) for k in range(d))

def get_all_union_sets(base_sets):
    """
    Generates all possible unique sets by taking unions of the base sets.
    A representation must have at least one character, so we exclude the empty set.
    """
    # The powerset recipe from the itertools documentation, starting from subsets of size 1.
    def powerset(iterable):
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))
    
    possible_unions = set()
    for subset_of_sets in powerset(base_sets):
        # Calculate the union of all sets in the current subset
        current_union = frozenset().union(*subset_of_sets)
        possible_unions.add(current_union)
    return possible_unions

# For the group G1 = Z_18
# The character orders are the divisors of 18.
divisors_18 = {1, 2, 3, 6, 9, 18}
char_images_g1 = {get_roots(d) for d in divisors_18}
union_sets_g1 = get_all_union_sets(char_images_g1)
num_sets_g1 = len(union_sets_g1)

# For the group G2 = Z_2 x Z_3 x Z_3
# The character orders are the divisors of the exponent, lcm(2, 3, 3) = 6.
divisors_6 = {1, 2, 3, 6}
char_images_g2 = {get_roots(d) for d in divisors_6}
union_sets_g2 = get_all_union_sets(char_images_g2)
num_sets_g2 = len(union_sets_g2)

# Find the number of common sets using set intersection
num_common_sets = len(union_sets_g1.intersection(union_sets_g2))

# Calculate the total number of unique sets using the principle of inclusion-exclusion
total_unique_sets = num_sets_g1 + num_sets_g2 - num_common_sets

# Print the breakdown of the final calculation
print(f"Number of unique sets for Z_18: {num_sets_g1}")
print(f"Number of unique sets for Z_2 x Z_3 x Z_3: {num_sets_g2}")
print(f"Number of common sets: {num_common_sets}")
print("Total number of unique sets is the size of the union of the two collections of sets.")
print(f"Final calculation: {num_sets_g1} + {num_sets_g2} - {num_common_sets} = {total_unique_sets}")
print(f"\nThus, there are {total_unique_sets} unique sets S(p) ∩ D.")
