def solve_alien_invasion():
    """
    Solves the alien invasion puzzle by simulating the growth of the colony
    from an optimized initial configuration.
    """
    # The board is represented as a set of occupied coordinates (col, row).
    # We use 0-indexed coordinates, where (0,0) is a1, (7,7) is h8.
    # 'd5' corresponds to (3, 4), 'e5' to (4, 4).

    # The optimal strategy found is a 3-line staggered configuration.
    # The 8 initial squares are distributed as follows:
    # Line 1 (row 2, i.e., the '3' rank): 1 square at d3 -> (3, 2)
    # Line 2 (row 4, i.e., the '5' rank): 2 squares at d5, e5 -> (3, 4), (4, 4)
    # Line 3 (row 6, i.e., the '7' rank): 5 squares at d7-h7 -> (3,6) to (7,6)
    initial_squares = {
        (3, 2),
        (3, 4), (4, 4),
        (3, 6), (4, 6), (5, 6), (6, 6), (7, 6)
    }

    occupied = set(initial_squares)
    
    turn_captures_counts = []
    
    while True:
        capturable = set()
        # Iterate over all squares on the 8x8 board
        for r in range(8):
            for c in range(8):
                if (c, r) in occupied:
                    continue

                neighbors = 0
                # Check 4 neighbors (up, down, left, right)
                if (c, r - 1) in occupied: neighbors += 1
                if (c, r + 1) in occupied: neighbors += 1
                if (c - 1, r) in occupied: neighbors += 1
                if (c + 1, r) in occupied: neighbors += 1
                
                if neighbors >= 2:
                    capturable.add((c, r))
        
        if not capturable:
            break
            
        turn_captures_counts.append(len(capturable))
        occupied.update(capturable)

    initial_size = len(initial_squares)
    final_size = len(occupied)
    
    equation_parts = [str(initial_size)] + [str(n) for n in turn_captures_counts]
    
    print("The maximal size K is found by simulating an optimal initial arrangement.")
    print("The process unfolds as follows:")
    print(f"Initial colony size: {initial_size}")
    for i, count in enumerate(turn_captures_counts):
        print(f"Squares captured in turn {i+1}: {count}")
    
    print("\nThe final size K is the sum of all squares:")
    print(f"K = {' + '.join(equation_parts)} = {final_size}")

solve_alien_invasion()
<<<25>>>