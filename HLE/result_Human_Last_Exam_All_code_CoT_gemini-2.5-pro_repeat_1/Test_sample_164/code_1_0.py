import sys

def solve_reversi():
    """
    Calculates the maximum number of disks Black can flip in a single turn
    on the given Reversi board.
    """
    # Board dimensions
    ROWS, COLS = 20, 20

    # Player and opponent symbols
    PLAYER = 'B'
    OPPONENT = 'W'
    EMPTY = '.'

    # Initialize an empty 20x20 board
    board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]

    # Helper function to convert alphanumeric coordinates (e.g., "A1", "T20")
    # to 0-indexed (row, col) tuples.
    def to_coords(s):
        row_char = s[0]
        col_num = int(s[1:])
        # 'A' becomes row 0, 'B' becomes row 1, etc.
        # Column 1 becomes col 0, column 2 becomes col 1, etc.
        return ord(row_char) - ord('A'), col_num - 1

    # List of initial disk positions from the problem description
    black_disks_pos = ["F9", "G9", "H9", "I9", "I13", "J8", "J9", "J15", "M9", "N9", "O10", "O11"]
    white_disks_pos = [
        "I10", "I11", "I12", "J10", "J11", "J12", "J13", "J14", "K8", "K9",
        "K10", "K11", "K12", "L8", "L9", "L10", "L11", "L12", "M10", "M11",
        "M12", "N10", "N11", "N12", "O12"
    ]

    # Populate the board with the disks
    for pos_str in black_disks_pos:
        r, c = to_coords(pos_str)
        board[r][c] = PLAYER

    for pos_str in white_disks_pos:
        r, c = to_coords(pos_str)
        board[r][c] = OPPONENT

    # The 8 directions (row_change, col_change) to check for flips
    directions = [(-1, -1), (-1, 0), (-1, 1),
                  (0, -1),          (0, 1),
                  (1, -1), (1, 0), (1, 1)]

    max_flips = 0

    # Iterate over every cell on the board to find all possible moves
    for r_start in range(ROWS):
        for c_start in range(COLS):
            # A move can only be made on an empty square
            if board[r_start][c_start] != EMPTY:
                continue

            flips_for_this_move = 0
            is_valid_move = False
            
            # For each empty square, check all 8 directions for potential flips
            for dr, dc in directions:
                r, c = r_start + dr, c_start + dc
                
                disks_in_line = []
                
                # Trace along the current direction
                while 0 <= r < ROWS and 0 <= c < COLS and board[r][c] == OPPONENT:
                    disks_in_line.append((r, c))
                    r += dr
                    c += dc

                # If the line of opponent disks is terminated by one of our own disks,
                # it's a valid flip sequence.
                if 0 <= r < ROWS and 0 <= c < COLS and board[r][c] == PLAYER:
                    if len(disks_in_line) > 0:
                        is_valid_move = True
                        flips_for_this_move += len(disks_in_line)

            # Update the maximum number of flips if this move is the best so far
            if is_valid_move and flips_for_this_move > max_flips:
                max_flips = flips_for_this_move

    # Print the final result
    print("The largest number of disks that black can flip is:")
    print(max_flips)

solve_reversi()