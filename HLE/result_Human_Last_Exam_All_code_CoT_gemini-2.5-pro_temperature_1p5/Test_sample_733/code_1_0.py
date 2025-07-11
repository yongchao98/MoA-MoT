import itertools

def check_groups():
    """
    This program finds the number of finite groups containing a maximal 
    product-free set of size 2 by computationally checking small groups.
    """

    # A helper class to represent a finite group using its Cayley table
    class Group:
        def __init__(self, name, elements, mul_table):
            self.name = name
            # Elements are represented by their indices 0 to n-1
            self.num_elements = len(elements)
            self.mul_table = mul_table

        def multiply(self, idx1, idx2):
            return self.mul_table[idx1][idx2]

        def has_mpfs2(self):
            """Checks if the group has a maximal product-free set of size 2."""
            element_indices = list(range(self.num_elements))
            # The identity element (index 0 by convention) cannot be in a product-free set
            non_identity_indices = element_indices[1:]

            # Iterate over all subsets S of size 2
            for s_indices in itertools.combinations(non_identity_indices, 2):
                s = set(s_indices)
                
                # 1. Check if S is product-free
                is_product_free = True
                product_set = {self.multiply(i, j) for i in s for j in s}
                if not s.isdisjoint(product_set):
                    is_product_free = False
                
                if not is_product_free:
                    continue

                # 2. If S is product-free, check if it's maximal
                is_maximal = True
                other_indices = [idx for idx in element_indices if idx not in s]
                
                for g_idx in other_indices:
                    t = s.union({g_idx})
                    
                    # If S U {g} is also product-free, then S is not maximal
                    t_product_set = {self.multiply(i, j) for i in t for j in t}
                    if t.isdisjoint(t_product_set):
                        is_maximal = False
                        break
                
                if is_maximal:
                    # Found one! This group qualifies.
                    return True
            
            return False

    # --- Group Definitions ---
    def make_cyclic(n):
        name = f"C_{n}"
        elements = list(range(n))
        table = [[(i + j) % n for j in elements] for i in elements]
        return Group(name, elements, table)

    def make_dihedral(n): # D_{2n}
        order = 2 * n
        name = f"D_{order}"
        elements = [(i, j) for i in range(n) for j in range(2)] # (rotation, reflection_flag)
        elem_to_idx = {el: i for i, el in enumerate(elements)}
        table = [[0]*order for _ in range(order)]
        for idx1, (i1, j1) in enumerate(elements):
            for idx2, (i2, j2) in enumerate(elements):
                if j1 == 0:
                    res_i, res_j = (i1 + i2) % n, j2
                else:
                    res_i, res_j = (i1 - i2 + n) % n, (1 + j2) % 2
                table[idx1][idx2] = elem_to_idx[(res_i, res_j)]
        return Group(name, elements, table)

    def make_direct_product(g1, g2):
        name = f"{g1.name} x {g2.name}"
        n1, n2 = g1.num_elements, g2.num_elements
        order = n1 * n2
        elements = [(i,j) for i in range(n1) for j in range(n2)]
        elem_to_idx = {el: i for i, el in enumerate(elements)}
        table = [[0]*order for _ in range(order)]
        for idx1, (i1, j1) in enumerate(elements):
            for idx2, (i2, j2) in enumerate(elements):
                res_i = g1.multiply(i1, i2)
                res_j = g2.multiply(j1, j2)
                table[idx1][idx2] = elem_to_idx[(res_i, res_j)]
        return Group(name, elements, table)

    # --- Main Logic ---
    groups_to_check = [
        # Check all groups up to order 10
        make_cyclic(2), make_cyclic(3), make_cyclic(4), make_cyclic(5), 
        make_cyclic(6), make_cyclic(7), make_cyclic(8), make_cyclic(9), make_cyclic(10),
        make_direct_product(make_cyclic(2), make_cyclic(2)), # Klein four-group
        make_direct_product(make_cyclic(2), make_cyclic(4)),
        make_direct_product(make_cyclic(3), make_cyclic(3)),
        make_dihedral(3), # D_6 (S_3)
        make_dihedral(4), # D_8
        make_dihedral(5), # D_10
    ]
    
    found_groups = []
    for g in groups_to_check:
        # Renaming C_2 x C_2 for clarity, as it's commonly known as V_4
        group_name = g.name
        if group_name == 'C_2 x C_2':
            group_name = 'C_2 x C_2 (Klein group V_4)'

        if g.has_mpfs2():
            found_groups.append(group_name)
    
    found_groups.sort()

    print("The finite groups known to contain a maximal product-free set of size 2 are:")
    for name in found_groups:
        print(f"- {name}")
    
    print("\nTo find the total number of these groups, we count them:")
    
    equation_parts = ["1" for _ in found_groups]
    equation_str = " + ".join(equation_parts)
    
    print(f"The count is: {equation_str} = {len(found_groups)}")

check_groups()