def solve():
    """
    Calculates the number of isomorphism classes of vertex-transitive graphs
    with 8 vertices for each degree j from 0 to 7.
    """
    # The total number of non-isomorphic vertex-transitive graphs on 8 vertices.
    total_graphs = 22

    # Initialize the list of counts n_j for j=0..7
    counts = [0] * 8

    # Step 4: Determine counts for small degrees by direct enumeration.
    # n_0: The empty graph (8 isolated vertices).
    counts[0] = 1
    # n_1: The perfect matching (4 disjoint edges).
    counts[1] = 1
    # n_2: The 8-cycle (C8) and the disjoint union of two 4-cycles (2C4).
    counts[2] = 2

    # Step 5: Use the complementation property (n_j = n_{7-j}) for high degrees.
    counts[7] = counts[0]
    counts[6] = counts[1]
    counts[5] = counts[2]

    # Step 6: Calculate the sum of known counts.
    known_sum = sum(counts)
    
    # Calculate the sum of remaining counts (for degrees 3 and 4).
    remaining_sum = total_graphs - known_sum

    # Step 7: The remaining counts n_3 and n_4 must be equal due to symmetry.
    counts[3] = remaining_sum // 2
    counts[4] = remaining_sum // 2
    
    # Print the final result in the specified format.
    print(counts)

solve()