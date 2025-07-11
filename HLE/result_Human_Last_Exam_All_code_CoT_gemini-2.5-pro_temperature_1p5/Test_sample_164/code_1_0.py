import collections

def parse_coords(coord_str):
    """Converts string coordinates like 'F9' to 0-indexed (row, col) tuples."""
    row_char = coord_str[0]
    col_num = int(coord_str[1:])
    row = ord(row_char) - ord('A')
    col = col_num - 1
    return row, col

def setup_board():
    """Initializes the 20x20 Reversi board with pieces from the problem."""
    rows, cols = 20, 20
    board = [['.' for _ in range(cols)] for _ in range(rows)]

    black_pieces_str = "F9 G9 H9 I9 I13 J8 J9 J15 M9 N9 O10 O11"
    white_pieces_str = "I10 I11 I12 J10 J11 J12 J13 J14 K8 K9 K10 K11 K12 L8 L9 L10 L11 L12 M10 M11 M12 N10 N11 N12 O12"

    for s in black_pieces_str.split():
        r, c = parse_coords(s)
        if 0 <= r < rows and 0 <= c < cols:
            board[r][c] = 'B'

    for s in white_pieces_str.split():
        r, c = parse_coords(s)
        if 0 <= r < rows and 0 <= c < cols:
            board[r][c] = 'W'
            
    return board

def calculate_flips_for_move(board, r, c, player):
    """
    Calculates the number of opponent's disks flipped for a move at (r, c).
    Returns the total flips and a list of flips from each successful direction.
    """
    if board[r][c] != '.':
        return 0, []

    rows, cols = len(board), len(board[0])
    opponent = 'W' if player == 'B' else 'B'
    
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    
    total_flips = 0
    directional_flips = []

    for dr, dc in directions:
        line_to_flip = []
        curr_r, curr_c = r + dr, c + dc

        # Move along the direction as long as we see opponent's disks
        while 0 <= curr_r < rows and 0 <= curr_c < cols and board[curr_r][curr_c] == opponent:
            line_to_flip.append((curr_r, curr_c))
            curr_r += dr
            curr_c += dc

        # If the line ends with the player's own disk, the flips are valid
        if 0 <= curr_r < rows and 0 <= curr_c < cols and board[curr_r][curr_c] == player:
            if line_to_flip:
                num_flips = len(line_to_flip)
                total_flips += num_flips
                directional_flips.append(num_flips)
    
    return total_flips, directional_flips

def solve_reversi():
    """
    Finds the best move for Black and prints the maximum number of flips.
    """
    board = setup_board()
    rows, cols = len(board), len(board[0])
    
    max_flips = 0
    best_move_details = []

    for r in range(rows):
        for c in range(cols):
            if board[r][c] == '.':
                flips, details = calculate_flips_for_move(board, r, c, 'B')
                if flips > max_flips:
                    max_flips = flips
                    best_move_details = details

    # Sort details for a consistent output equation
    best_move_details.sort(reverse=True)
    
    equation_parts = [str(x) for x in best_move_details]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {max_flips}")

solve_reversi()