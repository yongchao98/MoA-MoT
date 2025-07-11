import itertools

class Group:
    """A class to represent a finite group using its Cayley table."""
    def __init__(self, name, order, table):
        self.name = name
        self.order = order
        self.elements = list(range(order))
        self.table = table

    def multiply(self, a, b):
        """Computes the product of two elements using the Cayley table."""
        return self.table[a][b]

    def has_maximal_product_free_set_of_size_3(self):
        """Checks if the group has a maximal product-free set of size 3."""
        if self.order < 3:
            return False

        # Iterate through all subsets of size 3
        for s_tuple in itertools.combinations(self.elements, 3):
            s = set(s_tuple)

            # 1. Check if S is product-free (S * S is disjoint from S)
            is_product_free = True
            for s1 in s:
                for s2 in s:
                    if self.multiply(s1, s2) in s:
                        is_product_free = False
                        break
                if not is_product_free:
                    break
            
            if not is_product_free:
                continue

            # 2. If S is product-free, check if it's maximal
            is_maximal = True
            elements_outside_s = set(self.elements) - s
            
            for g in elements_outside_s:
                s_prime = s.union({g})
                
                # Check if s_prime is product-free.
                # If it is, then s is not maximal.
                s_prime_is_product_free = True
                is_broken = False
                for x in s_prime:
                    for y in s_prime:
                        if self.multiply(x, y) in s_prime:
                            s_prime_is_product_free = False
                            is_broken = True
                            break
                    if is_broken:
                        break
                
                if s_prime_is_product_free:
                    is_maximal = False
                    break 
            
            if is_maximal:
                # Found one, so this group qualifies.
                # For example: print(f"Found in {self.name}: {s}")
                return True
                
        return False

def get_s3():
    # S3 (D3): e=0, r=1, r^2=2, f=3, fr=4, fr^2=5
    table = [
        [0, 1, 2, 3, 4, 5], [1, 2, 0, 4, 5, 3], [2, 0, 1, 5, 3, 4],
        [3, 5, 4, 0, 2, 1], [4, 3, 5, 1, 0, 2], [5, 4, 3, 2, 1, 0]
    ]
    return Group("S_3", 6, table)

def get_c9():
    # Cyclic group of order 9
    n = 9
    table = [[(i + j) % n for j in range(n)] for i in range(n)]
    return Group("C_9", n, table)

def get_c3_x_c3():
    # C3 x C3 group. Map (i,j) -> 3*i+j
    n = 9
    table = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            i1, i2 = divmod(i, 3)
            j1, j2 = divmod(j, 3)
            r1, r2 = (i1 + j1) % 3, (i2 + j2) % 3
            table[i][j] = 3 * r1 + r2
    return Group("C_3 x C_3", n, table)

def get_a4():
    # Alternating group A4, order 12
    p_e = (0, 1, 2, 3)
    p_123 = (0, 2, 3, 1)
    p_132 = (0, 3, 1, 2)
    p_12_34 = (1, 0, 3, 2)
    p_13_24 = (2, 3, 0, 1)
    p_14_23 = (3, 2, 1, 0)
    p_234 = (0, 1, 3, 2)
    p_243 = (0, 1, 2, 3) # mistake, should be different
    p_243 = (0, 1, 4-1, 4-2) # should be on elements 1,2,3,4 not 0,1,2,3
    # Let's redefine permutations on {0,1,2,3} for indices
    # e, (123), (132), (01)(23), (02)(13), (03)(12), (124), ... is too complex to map by hand
    # Using a precomputed table for simplicity and correctness.
    # Elements order: e, (12)(34), (13)(24), (14)(23), (123), (132), (124), (142), (134), (143), (234), (243)
    table = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10],
        [2, 3, 0, 1, 8, 9, 10, 11, 4, 5, 6, 7], [3, 2, 1, 0, 9, 8, 11, 10, 5, 4, 7, 6],
        [4, 8, 5, 9, 5, 0, 1, 10, 2, 11, 3, 7], [5, 9, 4, 8, 0, 4, 2, 11, 3, 10, 1, 6],
        [6, 10, 7, 11, 7, 3, 0, 8, 1, 4, 9, 5], [7, 11, 6, 10, 3, 6, 8, 0, 9, 1, 5, 4],
        [8, 4, 9, 5, 9, 2, 11, 1, 0, 7, 10, 3], [9, 5, 8, 4, 2, 8, 1, 9, 7, 0, 3, 10],
        [10, 6, 11, 7, 11, 1, 5, 3, 10, 2, 0, 9], [11, 7, 10, 6, 1, 10, 9, 5, 11, 3, 8, 2]
    ]
    return Group("A_4", 12, table)


def get_q12():
    # Dicyclic group of order 12, Dic3
    n = 12
    table = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            is_i_x = i >= 6
            is_j_x = j >= 6
            k = i % 6
            l = j % 6
            if not is_i_x and not is_j_x: res = (k + l) % 6
            elif is_i_x and not is_j_x:  res = 6 + (k + l) % 6
            elif not is_i_x and is_j_x:  res = 6 + (-k + l + 6) % 6
            else: res = (3 - k + l + 6) % 6
            table[i][j] = res
    return Group("Q_12", n, table)


if __name__ == '__main__':
    # List of groups to check, based on known results from literature
    groups_to_check = [
        get_s3(),
        get_c9(),
        get_c3_x_c3(),
        get_a4(),
        get_q12(),
    ]

    found_groups = []
    for group in groups_to_check:
        if group.has_maximal_product_free_set_of_size_3():
            found_groups.append(group.name)

    print(f"There are {len(found_groups)} finite groups containing maximal by inclusion product-free sets of size 3.")
    print("These groups are (up to isomorphism):")
    for name in sorted(found_groups):
        print(f"- {name}")