import collections

def solve_reversi():
    """
    Calculates the largest number of disks black can flip in one turn on the given Reversi board.
    """
    # Board dimensions
    ROWS, COLS = 20, 20

    # Player representation
    EMPTY, BLACK, WHITE = 0, 1, -1

    # Initialize the board
    board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]

    # Helper function to convert board coordinates (e.g., 'A1', 'T20') to (row, col) indices
    def parse_coords(s):
        row_char = s[0]
        col_num = int(s[1:])
        row = ord(row_char) - ord('A')
        col = col_num - 1
        return row, col

    # List of initial piece positions
    black_pieces_coords = ['F9', 'G9', 'H9', 'I9', 'I13', 'J8', 'J9', 'J15', 'M9', 'N9', 'O10', 'O11']
    white_pieces_coords = ['I10', 'I11', 'I12', 'J10', 'J11', 'J12', 'J13', 'J14', 'K8', 'K9', 'K10', 'K11',
                           'K12', 'L8', 'L9', 'L10', 'L11', 'L12', 'M10', 'M11', 'M12', 'N10', 'N11', 'N12', 'O12']

    # Populate the board
    for s in black_pieces_coords:
        r, c = parse_coords(s)
        board[r][c] = BLACK
    for s in white_pieces_coords:
        r, c = parse_coords(s)
        board[r][c] = WHITE

    # All 8 directions for checking flips
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

    max_flips = 0

    # Iterate over every square on the board to find the best move
    for r in range(ROWS):
        for c in range(COLS):
            # We can only place a disk on an empty square
            if board[r][c] != EMPTY:
                continue

            total_flips_for_move = 0
            is_valid_move = False

            # Check all 8 directions from the current empty square
            for dr, dc in directions:
                flips_in_direction = []
                # Start from the adjacent square in the current direction
                next_r, next_c = r + dr, c + dc

                # Trace a line of opponent's pieces
                while 0 <= next_r < ROWS and 0 <= next_c < COLS:
                    if board[next_r][next_c] == WHITE:
                        # Found an opponent's piece, add it to the potential flip list
                        flips_in_direction.append((next_r, next_c))
                        next_r += dr
                        next_c += dc
                    elif board[next_r][next_c] == BLACK:
                        # Found our own piece, this is a valid capture line
                        if flips_in_direction:
                            total_flips_for_move += len(flips_in_direction)
                            is_valid_move = True
                        break  # Stop tracing in this direction
                    else:  # EMPTY or out of bounds
                        break # Stop tracing in this direction

            # Update the maximum flips found so far
            if is_valid_move and total_flips_for_move > max_flips:
                max_flips = total_flips_for_move

    print(f"The largest number of disks that black can flip in one turn is: {max_flips}")

solve_reversi()