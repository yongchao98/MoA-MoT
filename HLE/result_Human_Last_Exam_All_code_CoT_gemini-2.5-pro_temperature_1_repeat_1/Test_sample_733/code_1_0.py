import itertools

class Group:
    """A simple class to represent a finite group with a Cayley table."""
    def __init__(self, name, elements, op_table):
        self.name = name
        self.elements = elements
        self.op_table = op_table
        self.identity = self._find_identity()

    def _find_identity(self):
        """Find the identity element of the group."""
        for e in self.elements:
            is_identity = True
            for g in self.elements:
                if self.op_table[e][g] != g or self.op_table[g][e] != g:
                    is_identity = False
                    break
            if is_identity:
                return e
        return None

def check_group(group):
    """
    Checks if a group contains a maximal by inclusion product-free set of size 2.
    """
    n = len(group.elements)
    
    # Iterate over all possible 2-element subsets S = {a, b}
    for a, b in itertools.combinations(group.elements, 2):
        # The identity element cannot be in a product-free set
        if a == group.identity or b == group.identity:
            continue

        S = {a, b}

        # 1. Check if S is product-free
        # The product of any two elements in S cannot be in S
        products_from_S = {
            group.op_table[a][a],
            group.op_table[b][b],
            group.op_table[a][b],
            group.op_table[b][a]
        }
        if not products_from_S.isdisjoint(S):
            continue

        # 2. Check if S is maximal by inclusion
        # For any g not in S, the set S U {g} must NOT be product-free
        is_maximal = True
        elements_outside_S = [g for g in group.elements if g not in S]
        
        for g in elements_outside_S:
            S_union_g = list(S) + [g]
            
            # Check if S_union_g is product-free. If it is, S is not maximal.
            s_union_g_is_pf = True
            for x1 in S_union_g:
                for x2 in S_union_g:
                    prod = group.op_table[x1][x2]
                    if prod in S_union_g:
                        s_union_g_is_pf = False
                        break
                if not s_union_g_is_pf:
                    break
            
            if s_union_g_is_pf:
                is_maximal = False
                break
        
        if is_maximal:
            # We found a set S for this group, so the condition is met.
            # We print one such set found.
            print(f"Group {group.name} has a maximal product-free set of size 2, for example: {{{a}, {b}}}")
            return True
            
    return False

def get_groups_to_check():
    """Returns a list of all non-isomorphic groups of order <= 8."""
    groups = []
    
    # Order 2: Z_2
    groups.append(Group("Z_2", [0,1], [[(i+j)%2 for j in range(2)] for i in range(2)]))
    # Order 3: Z_3
    groups.append(Group("Z_3", [0,1,2], [[(i+j)%3 for j in range(3)] for i in range(3)]))
    # Order 4: Z_4
    groups.append(Group("Z_4", [0,1,2,3], [[(i+j)%4 for j in range(4)] for i in range(4)]))
    # Order 4: Z_2 x Z_2 (Klein four-group)
    groups.append(Group("Z_2 x Z_2", [0,1,2,3], [[i^j for j in range(4)] for i in range(4)]))
    # Order 5: Z_5
    groups.append(Group("Z_5", [0,1,2,3,4], [[(i+j)%5 for j in range(5)] for i in range(5)]))
    # Order 6: Z_6
    groups.append(Group("Z_6", [0,1,2,3,4,5], [[(i+j)%6 for j in range(6)] for i in range(6)]))
    # Order 6: S_3 (Symmetric group)
    s3_table = [
        [0,1,2,3,4,5], [1,0,4,5,2,3], [2,5,0,4,3,1],
        [3,4,5,0,1,2], [4,3,1,2,5,0], [5,2,3,1,0,4]
    ]
    groups.append(Group("S_3", list(range(6)), s3_table))
    # Order 7: Z_7
    groups.append(Group("Z_7", list(range(7)), [[(i+j)%7 for j in range(7)] for i in range(7)]))
    # Order 8: Z_8, Q_8, D_4, Z_2xZ_4, Z_2xZ_2xZ_2
    groups.append(Group("Z_8", list(range(8)), [[(i+j)%8 for j in range(8)] for i in range(8)]))
    q8_table = [ # Quaternion Group
        [0,1,2,3,4,5,6,7], [1,2,3,0,5,6,7,4], [2,3,0,1,6,7,4,5], [3,0,1,2,7,4,5,6],
        [4,5,6,7,2,1,0,3], [5,6,7,4,1,0,3,2], [6,7,4,5,0,3,2,1], [7,4,5,6,3,2,1,0]
    ]
    groups.append(Group("Q_8", list(range(8)), q8_table))
    d4_table = [ # Dihedral Group D_4
        [0,1,2,3,4,5,6,7], [1,2,3,0,7,4,5,6], [2,3,0,1,6,7,4,5], [3,0,1,2,5,6,7,4],
        [4,5,6,7,0,1,2,3], [5,4,7,6,1,0,3,2], [6,7,4,5,2,3,0,1], [7,6,5,4,3,2,1,0]
    ]
    groups.append(Group("D_4", list(range(8)), d4_table))
    z2z4_table = [[((i//4+j//4)%2)*4 + (i%4+j%4)%4 for j in range(8)] for i in range(8)]
    groups.append(Group("Z_2 x Z_4", list(range(8)), z2z4_table))
    z2z2z2_table = [[i^j for j in range(8)] for i in range(8)]
    groups.append(Group("Z_2 x Z_2 x Z_2", list(range(8)), z2z2z2_table))
    
    return groups

if __name__ == '__main__':
    groups_to_check = get_groups_to_check()
    count = 0
    print("Checking finite groups for maximal product-free sets of size 2...")
    for group in groups_to_check:
        if check_group(group):
            count += 1
    
    print("\n--- Summary ---")
    print(f"The number of finite groups (up to isomorphism of order <= 8) that contain a maximal product-free set of size 2 is:")
    print(count)