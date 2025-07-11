import collections

def solve_reversi():
    """
    Calculates the maximum number of disks Black can flip in a single move on the given Reversi board.
    """
    BOARD_SIZE = 20
    EMPTY, BLACK, WHITE = '.', 'B', 'W'

    # Helper to convert chess-like coordinates to 0-indexed (row, col)
    def to_idx(coord):
        row_char = coord[0]
        col_num = int(coord[1:])
        return ord(row_char) - ord('A'), col_num - 1

    # Helper to convert 0-indexed (row, col) back to chess-like coordinates
    def to_coord(r, c):
        return f"{chr(ord('A') + r)}{c + 1}"

    # Initialize the board
    board = [[EMPTY for _ in range(BOARD_SIZE)] for _ in range(BOARD_SIZE)]

    # Positions of the pieces
    black_pieces = ["F9", "G9", "H9", "I9", "I13", "J8", "J9", "J15", "M9", "N9", "O10", "O11"]
    white_pieces = ["I10", "I11", "I12", "J10", "J11", "J12", "J13", "J14", "K8", "K9", "K10", "K11", "K12",
                    "L8", "L9", "L10", "L11", "L12", "M10", "M11", "M12", "N10", "N11", "N12", "O12"]

    for coord in black_pieces:
        r, c = to_idx(coord)
        board[r][c] = BLACK
    for coord in white_pieces:
        r, c = to_idx(coord)
        board[r][c] = WHITE

    directions = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
    
    max_flips = 0
    best_move_info = {}

    # Iterate over every square on the board
    for r_start in range(BOARD_SIZE):
        for c_start in range(BOARD_SIZE):
            # We can only place a disk on an empty square
            if board[r_start][c_start] == EMPTY:
                flips_for_current_move = []
                
                # Check all 8 directions for possible flips
                for dr, dc in directions:
                    line_to_flip = []
                    r, c = r_start + dr, c_start + dc

                    while 0 <= r < BOARD_SIZE and 0 <= c < BOARD_SIZE:
                        if board[r][c] == WHITE:
                            line_to_flip.append((r, c))
                        elif board[r][c] == BLACK:
                            # If we find a black disk, the line is flanked
                            if line_to_flip:
                                flips_for_current_move.append(len(line_to_flip))
                            break
                        else: # EMPTY square
                            break
                        r += dr
                        c += dc
                
                total_flips = sum(flips_for_current_move)
                if total_flips > max_flips:
                    max_flips = total_flips
                    best_move_info = {
                        "move": (r_start, c_start),
                        "breakdown": sorted([f for f in flips_for_current_move if f > 0], reverse=True)
                    }

    if max_flips > 0:
        breakdown_str = " + ".join(map(str, best_move_info["breakdown"]))
        print(f"{breakdown_str} = {max_flips}")
    else:
        print("No legal moves found for Black.")

solve_reversi()