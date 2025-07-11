import collections
from itertools import combinations

def check_group_for_mpfs(elements, multiply_func, identity, group_name):
    """
    Checks if a given group contains a maximal product-free set of size 3.
    
    Args:
        elements: A list of the group's elements.
        multiply_func: A function that takes two elements and returns their product.
        identity: The identity element of the group.
        group_name: The name of the group for printing.

    Returns:
        True if an MPFS of size 3 is found, False otherwise.
    """
    non_identity_elements = [e for e in elements if e != identity]

    if len(non_identity_elements) < 3:
        return False

    # 1. Iterate through all subsets of size 3
    for s_tuple in combinations(non_identity_elements, 3):
        S = set(s_tuple)
        
        # 2. Check if product-free
        is_product_free = True
        for s1 in S:
            for s2 in S:
                if multiply_func(s1, s2) in S:
                    is_product_free = False
                    break
            if not is_product_free:
                break
        
        if not is_product_free:
            continue

        # 3. If product-free, check for maximality
        is_maximal = True
        G_minus_S = [g for g in elements if g not in S]
        for g in G_minus_S:
            S_union_g = S.union({g})
            
            # For S to be maximal, for every g not in S, S U {g} must NOT be product-free.
            # So, if we find even one g where S U {g} IS product-free, S is not maximal.
            is_superset_product_free = True
            
            # Check products involving g. Products from S are already known to not be in S.
            # If g*g, s*g, or g*s lands in S_union_g, it's not product-free.
            found_product_in_superset = False
            for s_member in S:
                if multiply_func(s_member, g) in S_union_g: found_product_in_superset = True; break
                if multiply_func(g, s_member) in S_union_g: found_product_in_superset = True; break
            if not found_product_in_superset:
                if multiply_func(g, g) in S_union_g: found_product_in_superset = True
            
            if not found_product_in_superset:
                 # This S_union_g is product-free, so S is not maximal.
                 is_maximal = False
                 break
        
        if is_maximal:
            # We found one, so this group has the property.
            return True
            
    return False

def define_groups():
    """Defines the 6 base groups."""
    groups = []
    
    # Z_6
    groups.append(("Z_6", list(range(6)), lambda a, b: (a + b) % 6, 0))
    
    # S_3
    p = list(combinations(range(3), 2))
    s3_elems = ['e', (1,2), (1,3), (2,3), (1,2,3), (1,3,2)] # String reps
    from sympy.combinatorics import Permutation
    # Using sympy to compute products easily, but storing as string representation
    s3_perms = {
        'e': Permutation.identity(3), (1,2): Permutation(1,2), (1,3): Permutation(1,3), (2,3): Permutation(2,3),
        (1,2,3): Permutation(1,2,3), (1,3,2): Permutation(1,3,2)
    }
    s3_perm_to_str = {v: k for k,v in s3_perms.items()}
    def s3_op(a, b):
        p1 = s3_perms[a]
        p2 = s3_perms[b]
        return s3_perm_to_str[p1*p2]
    groups.append(("S_3", s3_elems, s3_op, 'e'))

    # Z_9
    groups.append(("Z_9", list(range(9)), lambda a, b: (a + b) % 9, 0))

    # Z_3 x Z_3
    z3z3_elems = [(i,j) for i in range(3) for j in range(3)]
    def z3z3_op(a,b): return ((a[0]+b[0])%3, (a[1]+b[1])%3)
    groups.append(("Z_3xZ_3", z3z3_elems, z3z3_op, (0,0)))

    # D_10 (Dihedral group of order 10)
    d10_elems = [('r', i) for i in range(5)] + [('s', i) for i in range(5)]
    def d10_op(a, b):
        t1, i1 = a
        t2, i2 = b
        if t1 == 'r' and t2 == 'r': return ('r', (i1 + i2) % 5)
        if t1 == 'r' and t2 == 's': return ('s', (i1 + i2) % 5)
        if t1 == 's' and t2 == 'r': return ('s', (i1 - i2 + 5) % 5)
        if t1 == 's' and t2 == 's': return ('r', (i1 - i2 + 5) % 5)
    groups.append(("D_10", d10_elems, d10_op, ('r',0)))

    # D_14 (Dihedral group of order 14)
    d14_elems = [('r', i) for i in range(7)] + [('s', i) for i in range(7)]
    def d14_op(a, b):
        t1, i1 = a; t2, i2 = b
        if t1 == 'r' and t2 == 'r': return ('r', (i1 + i2) % 7)
        if t1 == 'r' and t2 == 's': return ('s', (i1 + i2) % 7)
        if t1 == 's' and t2 == 'r': return ('s', (i1 - i2 + 7) % 7)
        if t1 == 's' and t2 == 's': return ('r', (i1 - i2 + 7) % 7)
    groups.append(("D_14", d14_elems, d14_op, ('r',0)))
    
    return groups

# Main execution
base_groups = define_groups()
count_of_valid_groups = 0

for name, elements, op, identity in base_groups:
    if check_group_for_mpfs(elements, op, identity, name):
        count_of_valid_groups += 1

# Final Answer
print(f"The number of minimal finite groups containing a maximal by inclusion product-free set of size 3 is {count_of_valid_groups}.")