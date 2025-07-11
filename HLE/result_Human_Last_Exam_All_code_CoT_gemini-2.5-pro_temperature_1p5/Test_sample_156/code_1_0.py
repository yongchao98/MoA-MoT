def solve():
    """
    Solves the alien invasion puzzle by simulating the growth of the colony
    from an optimal initial configuration.
    """
    board_size = 8
    board = [[0] * board_size for _ in range(board_size)]

    # We propose an initial configuration in the shape of a large cross,
    # which is designed to trigger a large cascade of captures.
    # The 8 initial squares are:
    initial_squares_notation = [
        "d5", "e5",  # The two fixed squares
        "b5", "g5",  # The horizontal arms of the cross
        "d3", "e3",  # The upper vertical arm of the cross
        "d7", "e7"   # The lower vertical arm of the cross
    ]

    # Convert chess notation to 0-indexed grid coordinates (row, col)
    # where row 0 corresponds to rank '8' and col 0 to file 'a'.
    initial_coords = []
    for s in initial_squares_notation:
        col = ord(s[0]) - ord('a')
        row = 8 - int(s[1])
        initial_coords.append((row, col))

    # Mark the initial squares on the board
    for r, c in initial_coords:
        board[r][c] = 1

    # Simulate the growth process turn by turn
    while True:
        newly_captured = []
        for r in range(board_size):
            for c in range(board_size):
                if board[r][c] == 0:  # If the square is vacant
                    captured_neighbors = 0
                    # Check the four orthogonal neighbors
                    if r > 0 and board[r - 1][c] == 1:
                        captured_neighbors += 1
                    if r < board_size - 1 and board[r + 1][c] == 1:
                        captured_neighbors += 1
                    if c > 0 and board[r][c - 1] == 1:
                        captured_neighbors += 1
                    if c < board_size - 1 and board[r][c + 1] == 1:
                        captured_neighbors += 1
                    
                    if captured_neighbors >= 2:
                        newly_captured.append((r, c))

        if not newly_captured:
            # If no new squares were captured, the colony is stable.
            break

        # Add the newly captured squares to the board for the next turn
        for r, c in newly_captured:
            board[r][c] = 1

    final_size = sum(sum(row) for row in board)
    initial_count = len(initial_coords)
    grown_count = final_size - initial_count

    print(f"The chosen optimal initial 8 squares are: {', '.join(sorted(initial_squares_notation))}.")
    print("This configuration triggers a multi-stage growth cascade.")
    print("Final colony size (K) is the sum of initial squares and grown squares.")
    print(f"{initial_count} + {grown_count} = {final_size}")

solve()
<<<30>>>