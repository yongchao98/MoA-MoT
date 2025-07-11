import math
from itertools import combinations, product
from functools import reduce

def get_divisors(n):
    """Computes all divisors of an integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def get_character_image_indices(group_factors):
    """
    For a group G = Z_n1 x Z_n2 x ..., computes the set of indices 'd'
    such that U_d (d-th roots of unity) is an image of a character of G.
    """
    if not group_factors:
        return {1}
    
    divisor_lists = [get_divisors(n) for n in group_factors]
    
    lcm_values = set()
    for combo in product(*divisor_lists):
        if combo:
            result_lcm = reduce(lcm, combo)
            lcm_values.add(result_lcm)
            
    return lcm_values

def find_unique_union_sets(indices):
    """
    Finds the number of unique unions of sets U_d for d in indices.
    This is equivalent to finding the number of distinct antichains in the
    poset (indices, |), where | denotes divisibility.
    """
    
    # Generate all non-empty subsets of indices
    all_subsets = []
    for r in range(1, len(indices) + 1):
        all_subsets.extend(combinations(indices, r))
        
    # For each subset, find its maximal elements under divisibility
    # These maximal elements form an antichain that uniquely defines the union.
    unique_antichains = set()
    for subset in all_subsets:
        maximal_elements = []
        for i in subset:
            is_maximal = True
            for j in subset:
                if i != j and i % j == 0: # j divides i
                    is_maximal = False
                    break
            if is_maximal:
                maximal_elements.append(i)
        
        # Store the antichain as a frozenset to count unique ones
        unique_antichains.add(frozenset(sorted(maximal_elements)))
        
    return unique_antichains

def main():
    """
    Main function to solve the problem.
    """
    # The two non-isomorphic Abelian groups of order 18 are
    # G1 = Z_18 and G2 = Z_2 x Z_3 x Z_3
    groups = {
        "Z_18": [18],
        "Z_2 x Z_3 x Z_3": [2, 3, 3]
    }
    
    total_unique_sets = set()
    
    print("Analyzing Abelian groups of order 18:")
    for name, factors in groups.items():
        # Step 1: Find the set of indices 'd' for the character images U_d
        char_image_indices = get_character_image_indices(factors)
        
        # Step 2: Find all unique unions of these U_d sets
        union_sets = find_unique_union_sets(char_image_indices)
        
        print(f"\nFor group G = {name}:")
        print(f"  - The indices `d` of possible character images U_d are: {sorted(list(char_image_indices))}")
        print(f"  - Number of unique sets S(rho) = {len(union_sets)}")
        
        # Add the found sets (represented by antichains) to the total set
        total_unique_sets.update(union_sets)
        
    print("\n-----------------------------------------")
    print(f"Total number of unique sets S(rho) across all groups of order 18:")
    print(f"{len(total_unique_sets)}")
    
if __name__ == "__main__":
    main()
