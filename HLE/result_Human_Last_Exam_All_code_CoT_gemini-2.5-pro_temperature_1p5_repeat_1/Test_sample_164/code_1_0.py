def solve_reversi():
    """
    Analyzes a Reversi board to find the best move for the black player.
    The best move is defined as the one that flips the maximum number of white disks.
    """
    ROWS, COLS = 20, 20
    EMPTY, BLACK, WHITE = '.', 'B', 'W'

    # Initialize an empty 20x20 board
    board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]

    # Positions of the disks on the board (row, col) using 0-based indexing
    # A=0, B=1, ... and 1=0, 2=1, ...
    positions = {
        'B': [(5, 8), (6, 8), (7, 8), (8, 8), (8, 12), (9, 7), (9, 8), (9, 14),
              (12, 8), (13, 8), (14, 9), (14, 10)],
        'W': [(8, 9), (8, 10), (8, 11),
              (9, 9), (9, 10), (9, 11), (9, 12), (9, 13),
              (10, 7), (10, 8), (10, 9), (10, 10), (10, 11),
              (11, 7), (11, 8), (11, 9), (11, 10), (11, 11),
              (12, 9), (12, 10), (12, 11),
              (13, 9), (13, 10), (13, 11),
              (14, 11)]
    }

    for r, c in positions[BLACK]:
        board[r][c] = BLACK
    for r, c in positions[WHITE]:
        board[r][c] = WHITE

    player = BLACK
    opponent = WHITE
    
    # All 8 directions to check
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), 
                  (0, 1), (1, -1), (1, 0), (1, 1)]

    max_flips = 0
    best_move_info = None

    for r in range(ROWS):
        for c in range(COLS):
            # A move can only be made on an empty square
            if board[r][c] == EMPTY:
                current_flips_total = 0
                flips_per_direction = []

                # Check all 8 directions from the current empty square
                for dr, dc in directions:
                    line_to_flip = []
                    x, y = r + dr, c + dc
                    
                    # Trace a line in the current direction
                    while 0 <= x < ROWS and 0 <= y < COLS:
                        if board[x][y] == opponent:
                            line_to_flip.append((x, y))
                        elif board[x][y] == player:
                            # If we found an anchor piece, the move is valid in this direction
                            if line_to_flip:
                                current_flips_total += len(line_to_flip)
                                flips_per_direction.append(len(line_to_flip))
                            break # Stop searching in this direction
                        else: # Found an empty square or edge of board
                            break # Stop searching in this direction
                        x += dr
                        y += dc

                # If this move is the best so far, save its details
                if current_flips_total > max_flips:
                    max_flips = current_flips_total
                    best_move_info = {
                        "row": r,
                        "col": c,
                        "flips": flips_per_direction
                    }

    if best_move_info:
        # Build the equation string
        equation_parts = [str(f) for f in best_move_info["flips"]]
        equation = " + ".join(equation_parts) + f" = {max_flips}"
        print(f"The largest number of disks black can flip is {max_flips}.")
        
        # Convert 0-indexed coords to 1-indexed for the explanation
        move_row_label = chr(ord('A') + best_move_info['row'])
        move_col_label = best_move_info['col'] + 1
        print(f"This is achieved by placing a disk at {move_row_label}{move_col_label}.")
        print(f"The calculation is: {equation}")
    else:
        print("There are no valid moves for black.")

solve_reversi()