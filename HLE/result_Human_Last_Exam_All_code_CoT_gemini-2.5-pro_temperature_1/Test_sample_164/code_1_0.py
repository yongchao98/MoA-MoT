def solve_reversi():
    """
    Calculates the maximum number of disks Black can flip in a single move
    on the given Reversi board.
    """
    ROWS, COLS = 20, 20
    EMPTY, BLACK, WHITE = '.', 'B', 'W'

    # Initialize an empty 20x20 board
    board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]

    # Black pieces coordinates (row, col)
    # A=0, B=1, ... | 1=0, 2=1, ...
    black_pieces = [
        (5, 8), (6, 8), (7, 8), (8, 8), (8, 12), (9, 7), (9, 8),
        (9, 14), (12, 8), (13, 8), (14, 9), (14, 10)
    ]

    # White pieces coordinates (row, col)
    white_pieces = [
        (8, 9), (8, 10), (8, 11), (9, 9), (9, 10), (9, 11), (9, 12),
        (9, 13), (10, 7), (10, 8), (10, 9), (10, 10), (10, 11), (11, 7),
        (11, 8), (11, 9), (11, 10), (11, 11), (12, 9), (12, 10), (12, 11),
        (13, 9), (13, 10), (13, 11), (14, 11)
    ]

    for r, c in black_pieces:
        board[r][c] = BLACK
    for r, c in white_pieces:
        board[r][c] = WHITE

    def is_valid(r, c):
        return 0 <= r < ROWS and 0 <= c < COLS

    def calculate_flips(r_start, c_start):
        """Calculates flips for placing a black disk at (r_start, c_start)."""
        if board[r_start][c_start] != EMPTY:
            return 0

        total_flips = 0
        # Directions: N, NE, E, SE, S, SW, W, NW
        directions = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]

        for dr, dc in directions:
            r, c = r_start + dr, c_start + dc
            disks_to_flip = []

            # Move in the current direction
            while is_valid(r, c):
                if board[r][c] == WHITE:
                    disks_to_flip.append((r, c))
                elif board[r][c] == BLACK:
                    # Found a bracketing black piece, add the disks to the total
                    total_flips += len(disks_to_flip)
                    break
                else: # Empty square
                    break
                r, c = r + dr, c + dc
        
        return total_flips

    max_flips = 0
    
    # Iterate over every square on the board
    for r in range(ROWS):
        for c in range(COLS):
            # We can only place a disk on an empty square
            if board[r][c] == EMPTY:
                flips = calculate_flips(r, c)
                if flips > max_flips:
                    max_flips = flips
    
    # Let's break down the best move for clarity
    # Based on analysis, one of the best moves is at H11 (7, 10)
    # Let's calculate its components.
    # South: flips I11, J11, K11, L11, M11, N11 (6 disks) before hitting O11(B)
    flips_south = 6
    # Southwest: flips I10 (1 disk) before hitting J9(B)
    flips_southwest = 1
    
    # Another best move is at M13 (12, 12)
    # West: flips M12, M11, M10 (3 disks) before hitting M9(B)
    # Northwest: flips L12, K11, J10 (3 disks) before hitting I9(B)
    # Southwest: flips N12 (1 disk) before hitting O11(B)
    
    print("The largest number of disks that black can flip is the sum of flips from all valid directions for the best move.")
    print(f"For a move at H11, this is {flips_south} (South) + {flips_southwest} (Southwest).")
    print(f"{flips_south} + {flips_southwest} = {max_flips}")


solve_reversi()
<<<7>>>