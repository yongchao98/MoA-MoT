import collections

def solve_reversi():
    """
    Solves the Reversi problem by finding the best move for Black.
    """
    BOARD_SIZE = 20
    PLAYER = 'B'
    OPPONENT = 'W'
    EMPTY = '.'

    # Initialize the board
    board = [[EMPTY for _ in range(BOARD_SIZE)] for _ in range(BOARD_SIZE)]

    # Helper to place pieces using algebraic notation (e.g., 'A1')
    def place_piece(pos, piece):
        row_char = pos[0]
        col_num = int(pos[1:])
        row = ord(row_char) - ord('A')
        col = col_num - 1
        if 0 <= row < BOARD_SIZE and 0 <= col < BOARD_SIZE:
            board[row][col] = piece

    # Black pieces
    black_pieces = [
        "F9", "G9", "H9", "I9", "I13", "J8", "J9", "J15", "M9", "N9", "O10", "O11"
    ]
    # White pieces
    white_pieces = [
        "I10", "I11", "I12", "J10", "J11", "J12", "J13", "J14", "K9", "K10",
        "K11", "K12", "L9", "L10", "L11", "L12", "M10", "M11", "M12", "N10",
        "N11", "N12", "O12"
    ]

    for pos in black_pieces:
        place_piece(pos, PLAYER)
    for pos in white_pieces:
        place_piece(pos, OPPONENT)

    # All 8 directions (dr, dc)
    DIRECTIONS = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]

    max_flips = 0
    best_move_details = []

    # Iterate over every empty square on the board
    for r in range(BOARD_SIZE):
        for c in range(BOARD_SIZE):
            if board[r][c] == EMPTY:
                current_flips = 0
                current_move_details = []
                
                # Check all 8 directions from the current empty square
                for dr, dc in DIRECTIONS:
                    flips_in_direction = []
                    line = []
                    
                    # Move one step in the current direction
                    nr, nc = r + dr, c + dc
                    
                    # Keep moving as long as we see opponent pieces
                    while 0 <= nr < BOARD_SIZE and 0 <= nc < BOARD_SIZE and board[nr][nc] == OPPONENT:
                        line.append((nr, nc))
                        nr += dr
                        nc += dc
                    
                    # If the line ends with one of our pieces, the move is valid
                    if 0 <= nr < BOARD_SIZE and 0 <= nc < BOARD_SIZE and board[nr][nc] == PLAYER:
                        if line: # Ensure we actually flanked something
                            flips_in_direction.extend(line)

                # If this direction resulted in flips, record it
                if flips_in_direction:
                    num_flips_in_dir = len(flips_in_direction)
                    current_flips += num_flips_in_dir
                    current_move_details.append(num_flips_in_dir)
                
                # Check if this move is the new best move
                if current_flips > max_flips:
                    max_flips = current_flips
                    best_move_details = current_move_details

    # The problem asks to show the equation of the sum.
    # The details of the best move are stored in best_move_details.
    # The list might not be sorted, so sort for a consistent output.
    best_move_details.sort(reverse=True)
    
    # Format the output string as an equation
    equation_parts = [str(num) for num in best_move_details]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {max_flips}")


solve_reversi()
<<<7>>>