import collections

def solve_reversi():
    """
    Solves the Reversi problem for the given board state.
    Finds the move for black that flips the maximum number of white disks.
    """
    # Board dimensions
    ROWS, COLS = 20, 20
    EMPTY, BLACK, WHITE = '.', 'B', 'W'

    # Initialize an empty board
    board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]

    # Helper to convert 'A1' style coordinates to 0-indexed (row, col)
    def to_coords(s):
        row_char = s[0]
        col_num = int(s[1:])
        return ord(row_char) - ord('A'), col_num - 1

    # Positions of the pieces
    black_pieces = [
        "F9", "G9", "H9", "I9", "I13", "J8", "J9", "J15", "M9", "N9", "O10", "O11"
    ]
    white_pieces = [
        "I10", "I11", "I12", "J10", "J11", "J12", "J13", "J14", "K8", "K9",
        "K10", "K11", "K12", "L8", "L9", "L10", "L11", "L12", "M10", "M11",
        "M12", "N10", "N11", "N12", "O12"
    ]

    # Populate the board
    for pos_str in black_pieces:
        r, c = to_coords(pos_str)
        board[r][c] = BLACK
    for pos_str in white_pieces:
        r, c = to_coords(pos_str)
        board[r][c] = WHITE

    # Directions to check (8 directions)
    directions = [
        (-1, -1, "Up-Left"), (-1, 0, "Up"), (-1, 1, "Up-Right"),
        (0, -1, "Left"), (0, 1, "Right"),
        (1, -1, "Down-Left"), (1, 0, "Down"), (1, 1, "Down-Right")
    ]

    max_flips = 0
    best_move_info = None

    # Iterate through every square on the board
    for r in range(ROWS):
        for c in range(COLS):
            # A move can only be made on an empty square
            if board[r][c] == EMPTY:
                current_flips = 0
                flips_by_direction = collections.defaultdict(int)
                
                # Check all 8 directions from the current empty square
                for dr, dc, dir_name in directions:
                    line_to_flip = []
                    nr, nc = r + dr, c + dc

                    # Trace along the direction
                    while 0 <= nr < ROWS and 0 <= nc < COLS:
                        if board[nr][nc] == WHITE:
                            line_to_flip.append((nr, nc))
                        elif board[nr][nc] == BLACK:
                            # Found a bracketing black piece, the move is valid in this direction
                            if line_to_flip:
                                flips_in_dir = len(line_to_flip)
                                current_flips += flips_in_dir
                                flips_by_direction[dir_name] = flips_in_dir
                            break
                        else: # Empty square
                            break
                        nr, nc = nr + dr, c + dc
                
                # If this move is better than the best one found so far, update
                if current_flips > max_flips:
                    max_flips = current_flips
                    # Helper to convert 0-indexed coords back to 'A1' style
                    move_str = f"{chr(ord('A') + r)}{c + 1}"
                    best_move_info = {
                        "move": move_str,
                        "total_flips": current_flips,
                        "details": flips_by_direction
                    }

    # Output the results
    if best_move_info:
        print(f"The best move for black is at {best_move_info['move']}.")
        
        flip_counts = []
        for direction, count in best_move_info['details'].items():
            print(f"It flips {count} white disk(s) in the {direction} direction.")
            flip_counts.append(str(count))
        
        equation = " + ".join(flip_counts)
        print(f"\nThe final calculation for the total number of flips is:")
        print(f"{equation} = {best_move_info['total_flips']}")
        
        print("\nThe largest number of disks that black can flip in one turn is:")
        print(best_move_info['total_flips'])
    else:
        print("No valid moves found for black.")

solve_reversi()
<<<7>>>