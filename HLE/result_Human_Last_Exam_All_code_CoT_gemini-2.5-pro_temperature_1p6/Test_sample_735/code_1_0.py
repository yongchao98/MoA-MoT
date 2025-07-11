import itertools

# Helper class to represent a finite group
class FiniteGroup:
    """A simple class to represent a finite group using an operation table."""
    def __init__(self, name, elements, op_table, identity):
        self.name = name
        self.elements = elements
        self.op_table = op_table
        self.identity = identity
        self.elem_to_idx = {elem: i for i, elem in enumerate(elements)}
        self.idx_to_elem = {i: elem for i, elem in enumerate(elements)}

    def multiply(self, e1, e2):
        """Computes the product of two elements in the group."""
        i1 = self.elem_to_idx[e1]
        i2 = self.elem_to_idx[e2]
        return self.idx_to_elem[self.op_table[i1][i2]]

    def size(self):
        return len(self.elements)

def is_product_free(s_set, group):
    """Checks if a set S is product-free in group G."""
    for e1 in s_set:
        for e2 in s_set:
            product = group.multiply(e1, e2)
            if product in s_set:
                return False
    return True

def find_maximal_product_free_sets_of_size_3(group):
    """
    Finds if a group contains a maximal by inclusion product-free set of size 3.
    """
    if group.size() < 4:
        return False
        
    elements_without_identity = [e for e in group.elements if e != group.identity]
    if len(elements_without_identity) < 3:
        return False
    
    # Iterate over all 3-element subsets of G \ {e}
    for s_tuple in itertools.combinations(elements_without_identity, 3):
        s_set = set(s_tuple)
        
        # 1. Check if S is product-free
        if not is_product_free(s_set, group):
            continue

        # 2. Check if S is maximal by inclusion
        is_maximal = True
        other_elements = [e for e in group.elements if e not in s_set]
        
        for g in other_elements:
            s_prime = s_set.union({g})
            
            # If S' = S U {g} is product-free for any g, then S is not maximal.
            if is_product_free(s_prime, group):
                is_maximal = False
                break
        
        if is_maximal:
            # Found one! No need to check other subsets for this group.
            return True
            
    return False

# --- Group Definitions ---

def get_cyclic_group(n):
    """Creates the cyclic group Z_n of order n."""
    elements = list(range(n))
    op_table = [[(i + j) % n for j in range(n)] for i in range(n)]
    return FiniteGroup(f"Z_{n}", elements, op_table, 0)

def get_dihedral_group(n):
    """Creates the dihedral group D_2n of order 2n."""
    elements = list(range(2 * n))
    op_table = [[0] * (2 * n) for _ in range(2 * n)]
    # Elements 0 to n-1 are rotations r^i, n to 2n-1 are reflections sr^i
    for i in range(2 * n):
        for j in range(2 * n):
            if i < n and j < n:  # rot * rot = r^i * r^j = r^(i+j)
                op_table[i][j] = (i + j) % n
            elif i < n and j >= n: # rot * ref = r^i * sr^j = sr^(j-i)
                op_table[i][j] = (((j - n) - i + n) % n) + n
            elif i >= n and j < n: # ref * rot = sr^i * r^j = sr^(i+j)
                op_table[i][j] = (((i - n) + j) % n) + n
            elif i >= n and j >= n: # ref * ref = sr^i * sr^j = r^(j-i)
                op_table[i][j] = (((j - n) - (i - n) + n) % n)

    name = f"D_{2*n}"
    if 2*n == 6:
        # S_3 is isomorphic to D_6
        name = "S_3 (isomorphic to D_6)"
    return FiniteGroup(name, elements, op_table, 0)

def get_quaternion_group():
    """Creates the quaternion group Q_8."""
    elements = ['1', '-1', 'i', '-i', 'j', '-j', 'k', '-k']
    elem_map = {e: i for i, e in enumerate(elements)}
    
    products = {
        ('i', 'i'): '-1', ('i', 'j'): 'k', ('i', 'k'): '-j',
        ('j', 'i'): '-k', ('j', 'j'): '-1', ('j', 'k'): 'i',
        ('k', 'i'): 'j', ('k', 'j'): '-i', ('k', 'k'): '-1'
    }
    op_table = [[0] * 8 for _ in range(8)]
    for e1_str in elements:
        for e2_str in elements:
            e1_sign = -1 if e1_str.startswith('-') else 1
            e2_sign = -1 if e2_str.startswith('-') else 1
            e1_base = e1_str[-1]
            e2_base = e2_str[-1]

            res_str = ''
            if e1_base == '1': res_str = e2_str
            elif e2_base == '1': res_str = e1_str
            else:
                if e1_base == e2_base:
                    res_str = '-1' if e1_sign * e2_sign == 1 else '1'
                else:
                    # e.g., i*j=k
                    p_res = products[(e1_base, e2_base)]
                    p_sign = -1 if p_res.startswith('-') else 1
                    res_base = p_res[-1]
                    # Total sign from inputs and cross product rule (e.g. j*i = -k)
                    final_sign = e1_sign * e2_sign * p_sign
                    res_str = res_base if final_sign == 1 else '-' + res_base

            op_table[elem_map[e1_str]][elem_map[e2_str]] = elem_map[res_str]

    return FiniteGroup("Q_8", elements, op_table, '1')

def solve():
    """
    This script searches for finite groups with maximal product-free sets of size 3.
    It checks a catalog of small groups against the standard mathematical definitions.
    """
    
    groups_to_check = []
    # Cyclic groups Z_n for n from 4 to 12
    for n in range(4, 13):
        groups_to_check.append(get_cyclic_group(n))
        
    # Dihedral groups D_2n for n from 3 to 6 (orders 6, 8, 10, 12)
    for n in range(3, 7):
        groups_to_check.append(get_dihedral_group(n))
        
    # Other common non-abelian groups
    groups_to_check.append(get_quaternion_group())
    
    found_groups = []
    print("Searching for groups with maximal product-free sets of size 3...")
    for group in groups_to_check:
        if find_maximal_product_free_sets_of_size_3(group):
            found_groups.append(group.name)
            
    print("\nBased on the computational search, the following groups were found:")
    for name in found_groups:
        print(f"- {name}")
    print("\nThe total number of such groups found is:")
    print(len(found_groups))

solve()
