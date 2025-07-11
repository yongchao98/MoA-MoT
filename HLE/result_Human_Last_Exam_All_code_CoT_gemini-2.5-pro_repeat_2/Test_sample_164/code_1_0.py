import collections

def solve():
    """
    Calculates the largest number of disks Black can flip in one turn on a custom Reversi board.
    """
    WIDTH, HEIGHT = 20, 20
    EMPTY, BLACK, WHITE = '.', 'B', 'W'
    PLAYER = BLACK
    OPPONENT = WHITE

    # Initialize the board
    board = [[EMPTY for _ in range(WIDTH)] for _ in range(HEIGHT)]

    # Black's disks
    black_disks_coords = [
        "F9", "G9", "H9", "I9", "I13", "J8", "J9", "J15",
        "M9", "N9", "O10", "O11"
    ]
    # White's disks
    white_disks_coords = [
        "I10", "I11", "I12", "J10", "J11", "J12", "J13", "J14",
        "K8", "K9", "K10", "K11", "K12", "L8", "L9", "L10",
        "L11", "L12", "M10", "M11", "M12", "N10", "N11", "N12", "O12"
    ]

    def parse_coords(s):
        row_char = s[0]
        col_str = s[1:]
        row = ord(row_char) - ord('A')
        col = int(col_str) - 1
        return row, col

    for coord in black_disks_coords:
        r, c = parse_coords(coord)
        if 0 <= r < HEIGHT and 0 <= c < WIDTH:
            board[r][c] = BLACK

    for coord in white_disks_coords:
        r, c = parse_coords(coord)
        if 0 <= r < HEIGHT and 0 <= c < WIDTH:
            board[r][c] = WHITE

    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    max_flips = 0
    
    # Iterate over every empty square to check for valid moves
    for r_start in range(HEIGHT):
        for c_start in range(WIDTH):
            if board[r_start][c_start] == EMPTY:
                current_move_flips = 0
                # Check all 8 directions from the empty square
                for dr, dc in directions:
                    disks_to_flip = []
                    r, c = r_start + dr, c_start + dc
                    
                    # Move along the direction as long as we see opponent's disks
                    while 0 <= r < HEIGHT and 0 <= c < WIDTH and board[r][c] == OPPONENT:
                        disks_to_flip.append((r, c))
                        r, c = r + dr, c + dc
                    
                    # If the line ends with the player's disk, the flips are valid
                    if 0 <= r < HEIGHT and 0 <= c < WIDTH and board[r][c] == PLAYER:
                        current_move_flips += len(disks_to_flip)
                
                if current_move_flips > max_flips:
                    max_flips = current_move_flips

    print(f"The largest number of disks that black can flip in one turn is: {max_flips}")

solve()