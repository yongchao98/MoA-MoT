import math
from fractions import Fraction
from itertools import chain, combinations

def get_roots_of_unity(k):
    """
    Returns a frozenset representing the k-th roots of unity.
    We represent the root e^(2*pi*i*j/k) by the fraction j/k to avoid floating point issues.
    The set {Fraction(j, k) for j in range(k)} uniquely identifies U_k.
    """
    if k == 0:
        return frozenset()
    # Normalize fractions to their simplest form to correctly handle unions,
    # e.g., U_2 = {0/2, 1/2} -> {0/1, 1/2} and U_4 = {0/4, 1/4, 2/4, 3/4} -> {0/1, 1/4, 1/2, 3/4}.
    # The union will then correctly recognize that 1/2 is a shared element.
    # The Fraction class handles this normalization automatically.
    return frozenset(Fraction(j, k) for j in range(k))

def get_all_possible_unions(image_sets):
    """
    Given a list of sets, computes all possible non-empty unions of these sets.
    Returns a set containing the unique union sets (represented as frozensets).
    """
    unique_unions = set()
    # Generate all non-empty subsets of image_sets
    power_set = chain.from_iterable(combinations(image_sets, r) for r in range(1, len(image_sets) + 1))
    
    for subset in power_set:
        # Compute the union of all sets in the current subset
        current_union = frozenset().union(*subset)
        unique_unions.add(current_union)
        
    return unique_unions

# --- Analysis for G1 = Z_18 ---
# The character image sets are U_d for d dividing 18.
divisors_18 = [1, 2, 3, 6, 9, 18]
char_images_G1 = [get_roots_of_unity(d) for d in divisors_18]
s_rho_sets_G1 = get_all_possible_unions(char_images_G1)
num_sets_G1 = len(s_rho_sets_G1)

# --- Analysis for G2 = Z_2 x Z_3 x Z_3 ---
# The character image sets are U_d for d dividing the exponent, lcm(2,3,3)=6.
divisors_6 = [1, 2, 3, 6]
char_images_G2 = [get_roots_of_unity(d) for d in divisors_6]
s_rho_sets_G2 = get_all_possible_unions(char_images_G2)
num_sets_G2 = len(s_rho_sets_G2)

# --- Combine results ---
# Find the number of sets that are possible for both types of groups.
num_common_sets = len(s_rho_sets_G1.intersection(s_rho_sets_G2))

# The total number of unique sets is |S1 U S2| = |S1| + |S2| - |S1 intersect S2|
total_unique_sets = num_sets_G1 + num_sets_G2 - num_common_sets

print(f"For the group G1 = Z_18, there are {num_sets_G1} possible unique sets S(rho).")
print(f"For the group G2 = Z_2 x Z_3 x Z_3, there are {num_sets_G2} possible unique sets S(rho).")
print(f"The number of unique sets common to both groups is {num_common_sets}.")
print(f"The total number of unique sets across all Abelian groups of order 18 is calculated by the Inclusion-Exclusion Principle:")
print(f"Total = (Sets from G1) + (Sets from G2) - (Common Sets)")
print(f"Total = {num_sets_G1} + {num_sets_G2} - {num_common_sets} = {total_unique_sets}")
