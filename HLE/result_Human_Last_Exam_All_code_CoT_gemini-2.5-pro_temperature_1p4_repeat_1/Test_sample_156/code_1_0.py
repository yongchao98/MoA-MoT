def solve_alien_colony():
    """
    Simulates the alien colony growth on a chessboard to find the maximal size.
    """
    # Board dimensions
    BOARD_SIZE = 8

    # Initial captured squares. Using (row, col) from (0,0) to (7,7).
    # d5 -> row 3, col 3
    # e5 -> row 3, col 4
    # We hypothesize a configuration that can fill a 6x6 square.
    # The 6x6 square is from row 1 to 6, and col 1 to 6.
    # Initial squares are d5, e5, and 6 others placed on the boundary of this 6x6 square.
    initial_squares = {
        # Fixed squares
        (3, 3), (3, 4),  # d5, e5

        # Corners of the 6x6 bounding box (b7, g7, b2, g2 in chess notation)
        (1, 1), (1, 6),
        (6, 1), (6, 6),

        # Midpoints of the vertical sides of the 6x6 box (b4, g4)
        (4, 1), (4, 6)
    }

    captured = set(initial_squares)
    
    while True:
        newly_captured = set()
        # Iterate over all squares on the board
        for r in range(BOARD_SIZE):
            for c in range(BOARD_SIZE):
                if (r, c) not in captured:
                    # Count captured neighbors for a vacant square
                    neighbor_count = 0
                    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        nr, nc = r + dr, c + dc
                        if (nr, nc) in captured:
                            neighbor_count += 1
                    
                    # Check if the square can be captured
                    if neighbor_count >= 2:
                        newly_captured.add((r, c))

        if not newly_captured:
            # No new squares were captured, the colony stops growing
            break
        
        # Add the new squares to the colony
        captured.update(newly_captured)

    # The result is the total number of captured squares
    k = len(captured)
    
    # Print the result as requested.
    print(k)

solve_alien_colony()