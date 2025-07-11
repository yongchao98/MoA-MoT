def solve_reversal_puzzle():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    total_moves = 0
    num_elements = 100
    
    # We need to perform swaps (i, 99-i) for i from 0 to 49.
    # There are 50 such pairs.
    num_pairs = num_elements // 2
    
    moves_per_pair = []

    for i in range(num_pairs):
        p1 = i
        p2 = num_elements - 1 - i
        
        # Physical distance between the elements to be swapped
        distance = p2 - p1
        
        # The cost of covering a distance 'd' is min(d%5, 5 - d%5) because
        # we can use free +/- 5 steps. This gives the minimum number of +/- 1 steps.
        # This is the shortest path distance 'dist' in our position graph.
        dist = min(distance % 5, 5 - (distance % 5))
        
        # The number of moves to perform a swap between two positions with a shortest
        # path distance of 'dist' is 2*dist - 1 (for dist > 0).
        if dist == 0:
            moves = 0
        else:
            moves = 2 * dist - 1
            
        moves_per_pair.append(moves)
        total_moves += moves

    # Print the full equation showing the sum of moves for each pair
    equation_str = " + ".join(map(str, moves_per_pair))
    print("The minimum number of moves is the sum of moves for each of the 50 pairs:")
    print(f"{equation_str} = {total_moves}")

solve_reversal_puzzle()