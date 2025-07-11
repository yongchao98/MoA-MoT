import sympy
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.groups import PermutationGroup
from itertools import combinations

def solve():
    """
    This function verifies that the dihedral group D_10 is a filled group.
    A group G is filled if every maximal product-free set S in G generates G.
    A set S is product-free if for any s1, s2 in S, s1*s2 is not in S.
    A product-free set is maximal if it's not a proper subset of any other product-free set.
    """
    
    # Define parameters for the group of order 2*q^m
    q = 5
    m = 1
    order = 2 * (q**m)

    # Step 1: Define the group G = D_10 (D_{2*5}) using permutations
    # Generators for D_10: rotation 'a' and reflection 'b'
    a = Permutation(list(range(q)))  # Rotation (0 1 2 3 4)
    b = Permutation(1, 4)(2, 3)      # Reflection
    G = PermutationGroup([a, b])
    
    elements = list(G.generate())
    identity = Permutation(q-1)

    print(f"Analyzing the group G = D_{2*q^m} = D_{order}")
    print(f"The equation for the group order is: 2 * {q}^{m} = {order}")
    print(f"Group has {G.order()} elements.")

    # Function to check if a set is product-free
    def is_product_free(s_tuple):
        s_set = set(s_tuple)
        if identity in s_set:
            return False
        for p1 in s_tuple:
            for p2 in s_tuple:
                if (p1 * p2) in s_set:
                    return False
        return True

    # Step 2: Find all product-free sets. This is computationally intensive.
    # We generate subsets and test them.
    # Using frozenset to store unique sets of elements.
    product_free_sets = set()
    # Iterate through all possible subset sizes
    for i in range(1, len(elements) + 1):
        # Iterate through all subsets of that size
        for subset_tuple in combinations(elements, i):
            if is_product_free(subset_tuple):
                product_free_sets.add(frozenset(subset_tuple))
    
    print(f"Found {len(product_free_sets)} product-free sets.")

    # Step 3: Find the maximal product-free sets from the list of all product-free sets.
    maximal_product_free_sets = []
    product_free_list = sorted(list(product_free_sets), key=len)
    for i, s1 in enumerate(product_free_list):
        is_maximal = True
        # Check if s1 is a proper subset of any other product-free set s2.
        # We only need to check sets s2 that are larger than s1.
        for j in range(i + 1, len(product_free_list)):
            s2 = product_free_list[j]
            if len(s1) < len(s2) and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_product_free_sets.append(list(s1))
            
    print(f"Found {len(maximal_product_free_sets)} maximal product-free sets.")
    print("-" * 20)

    # Step 4: For each maximal set, check if it generates the group.
    all_generate = True
    for i, s in enumerate(maximal_product_free_sets):
        subgroup = PermutationGroup(s)
        generates_G = (subgroup.order() == G.order())
        
        # To avoid overly long output, we'll just show the size of the set
        # and whether it generates the full group.
        # print(f"Maximal product-free set {i+1} (size {len(s)}): {s}")
        print(f"Maximal product-free set #{i+1} (size {len(s)}) generates G? {generates_G}")

        if not generates_G:
            all_generate = False
            # We can stop early if we find one that doesn't generate G.
            break
            
    print("-" * 20)

    # Final conclusion based on the exhaustive check
    if all_generate:
        print(f"Conclusion: The group D_{order} is FILLED, as all its maximal product-free sets generate it.")
    else:
        print(f"Conclusion: The group D_{order} is NOT FILLED, as a maximal product-free set that does not generate it was found.")

solve()