import itertools
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.groups import PermutationGroup

def find_maximal_product_free_sets(group):
    """
    Finds all maximal by inclusion product-free sets for a given group.
    
    A set S is product-free if for all x, y in S, x*y is not in S.
    A product-free set S is maximal if it's not a proper subset of another
    product-free set.
    """
    elements = list(group.elements)
    identity = group.identity
    
    # We only need to consider subsets of non-identity elements
    non_identity_elements = [el for el in elements if el != identity]
    
    # Generate all non-empty subsets of the non-identity elements
    powerset = []
    for i in range(1, len(non_identity_elements) + 1):
        for subset_tuple in itertools.combinations(non_identity_elements, i):
            powerset.append(frozenset(subset_tuple))

    # Filter for product-free sets
    product_free_sets = []
    for s in powerset:
        is_pf = True
        for x in s:
            for y in s:
                if (x * y) in s:
                    is_pf = False
                    break
            if not is_pf:
                break
        if is_pf:
            product_free_sets.append(s)
            
    # Filter for maximal product-free sets
    maximal_sets = []
    for s1 in product_free_sets:
        is_maximal = True
        for s2 in product_free_sets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_sets.append(s1)
            
    return maximal_sets

def analyze_group(name, group):
    """
    Analyzes a group for the filled and nilpotent properties.
    """
    print(f"--- Analyzing Group: {name} ---")
    
    # 1. Check if the group is nilpotent
    is_nilpotent = group.is_nilpotent
    print(f"Is {name} nilpotent? {is_nilpotent}")
    
    # 2. Check if the group is 'filled'
    # (i.e., union of maximal product-free sets equals G - {e})
    elements = set(group.elements)
    identity = group.identity
    non_identity_elements = elements - {identity}
    
    mfs_list = find_maximal_product_free_sets(group)
    
    # Print the maximal product-free sets found
    print(f"Maximal Product-Free Sets in {name}:")
    if not mfs_list:
        print("  None")
    else:
        for s in mfs_list:
            # Print sets in a readable format
            print(f"  {{ {', '.join(str(p) for p in sorted(list(s)))} }}")

    union_of_mfs = set().union(*mfs_list)
    
    is_filled = (union_of_mfs == non_identity_elements)
    print(f"Is {name} filled (i.e., Union(MFS) == G-{{e}})? {is_filled}")
    
    # Return if it satisfies both conditions
    return is_filled and is_nilpotent

def main():
    """
    Main function to identify finite filled nilpotent groups.
    """
    # Define the groups to check based on the theorem by Das & Jana
    # C_2, cyclic group of order 2
    c2 = PermutationGroup(Permutation(1)) 
    # C_3, cyclic group of order 3
    c3 = PermutationGroup(Permutation(2))
    # S_3, symmetric group of order 6
    s3 = PermutationGroup(Permutation(2), Permutation(0, 1))

    groups_to_check = {
        "C2": c2,
        "C3": c3,
        "S3": s3
    }
    
    filled_nilpotent_groups = []
    
    for name, group in groups_to_check.items():
        if analyze_group(name, group):
            filled_nilpotent_groups.append(name)
        print("-" * (len(name) + 21) + "\n")
        
    print("==================================================")
    print("Conclusion:")
    if not filled_nilpotent_groups:
        print("There are no finite filled nilpotent groups among the candidates.")
    else:
        final_answer = ", ".join(sorted(filled_nilpotent_groups))
        print(f"The finite filled nilpotent groups are: {final_answer}")
        print("Each of these groups satisfies two conditions:")
        print("1. Nilpotent: The group is the direct product of its Sylow subgroups.")
        print("2. Filled: The union of its maximal product-free sets equals the set of all its non-identity elements.")


if __name__ == "__main__":
    main()
<<<C2, C3>>>