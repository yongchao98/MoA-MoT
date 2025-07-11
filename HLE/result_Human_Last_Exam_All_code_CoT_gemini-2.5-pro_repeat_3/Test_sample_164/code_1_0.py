import collections

def solve_reversi():
    """
    Calculates the maximum number of disks Black can flip in a single move
    on the given Reversi board.
    """
    board_str = [
        "....................",  # A
        "....................",  # B
        "....................",  # C
        "....................",  # D
        "....................",  # E
        "........B...........",  # F
        "........B...........",  # G
        "........B...........",  # H
        "........BWWWB.......",  # I
        ".......BBWWWWWB.....",  # J
        ".......WWWWW........",  # K
        ".......WWWWW........",  # L
        "........BWWW........",  # M
        "........BWWW........",  # N
        ".........BBW........",  # O
        "....................",  # P
        "....................",  # Q
        "....................",  # R
        "....................",  # S
        "....................",  # T
    ]

    board = [list(row) for row in board_str]
    rows, cols = len(board), len(board[0])
    player = 'B'
    opponent = 'W'

    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

    max_flips = 0

    for r in range(rows):
        for c in range(cols):
            # Check only empty squares for potential moves
            if board[r][c] == '.':
                current_flips = 0
                is_valid_move = False

                # Check all 8 directions from the current empty square
                for dr, dc in directions:
                    flips_in_line = []
                    # Move one step in the current direction
                    nr, nc = r + dr, c + dc

                    # Trace along the line
                    while 0 <= nr < rows and 0 <= nc < cols:
                        if board[nr][nc] == opponent:
                            flips_in_line.append((nr, nc))
                        elif board[nr][nc] == player:
                            # Found an anchor piece, the move is valid in this direction
                            if flips_in_line:
                                current_flips += len(flips_in_line)
                                is_valid_move = True
                            break
                        else: # Found an empty square or went off board
                            break
                        
                        nr, nc = nr + dr, nc + dc
                
                if is_valid_move:
                    if current_flips > max_flips:
                        max_flips = current_flips

    print("The largest number of disks that black can flip in one turn is:")
    print(max_flips)

solve_reversi()