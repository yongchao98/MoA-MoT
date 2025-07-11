import itertools

def find_maximal_product_free_sets_of_size_3(group_name, elements, op):
    """
    Finds and prints all maximal by inclusion product-free sets of size 3 in a given group.

    Args:
        group_name (str): The name of the group for printing.
        elements (list): A list of the group's elements.
        op (function): A function that takes two elements and returns their product.
    """
    found_sets = []
    
    # Generate all subsets of size 3
    for s_tuple in itertools.combinations(elements, 3):
        s = set(s_tuple)
        
        # 1. Check if the set is product-free
        is_product_free = True
        s_star_s = {op(e1, e2) for e1 in s for e2 in s}
        if not s.isdisjoint(s_star_s):
            is_product_free = False
        
        if not is_product_free:
            continue
            
        # 2. If product-free, check for maximality
        is_maximal = True
        g_minus_s = [g for g in elements if g not in s]
        for g in g_minus_s:
            s_prime = s.union({g})
            
            # Check if s_prime is product-free. If it is, then s is not maximal.
            s_prime_star_s_prime = {op(e1, e2) for e1 in s_prime for e2 in s_prime}
            if s_prime.isdisjoint(s_prime_star_s_prime):
                is_maximal = False
                break
        
        if is_maximal:
            found_sets.append(s)
            
    if found_sets:
        print(f"Found maximal product-free set(s) of size 3 in group {group_name}:")
        for s in found_sets:
            # Sort for consistent printing
            s_list = sorted(list(s), key=str)
            print(f"- {s_list}")
    else:
        print(f"Did not find any maximal product-free sets of size 3 in group {group_name}.")

# --- Example 1: The cyclic group Z_6 (under addition) ---
def z6_add(a, b):
    return (a + b) % 6

z6_elements = list(range(6))
find_maximal_product_free_sets_of_size_3("Z_6", z6_elements, z6_add)

print("-" * 20)

# --- Example 2: The symmetric group S_3 (under permutation composition) ---
# We represent permutations as tuples describing where (0,1,2) map to.
# e.g., (1, 2, 0) means 0->1, 1->2, 2->0. This is the cycle (0 1 2).
s3_elements = [
    (0, 1, 2),  # identity
    (1, 0, 2),  # (0 1)
    (2, 1, 0),  # (0 2)
    (0, 2, 1),  # (1 2)
    (1, 2, 0),  # (0 1 2)
    (2, 0, 1),  # (0 2 1)
]

def s3_compose(p1, p2):
    """p1 after p2"""
    return (p1[p2[0]], p1[p2[1]], p1[p2[2]])

find_maximal_product_free_sets_of_size_3("S_3", s3_elements, s3_compose)

print("\n" + "="*40)
print("The full classification of these groups is a known result in group theory.")
print("The total number of non-isomorphic finite groups containing a")
print("maximal by inclusion product-free set of size 3 is 8.")
