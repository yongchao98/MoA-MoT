import collections

def solve_reversi():
    """
    Calculates the largest number of disks Black can flip in one turn on the given Reversi board.
    """
    WIDTH, HEIGHT = 20, 20
    PLAYER, OPPONENT = 'B', 'W'

    # Initialize a 20x20 board with empty squares
    board = [['.' for _ in range(WIDTH)] for _ in range(HEIGHT)]

    # Helper to place disks using algebraic notation (e.g., 'A', 1)
    def set_disk(row_char, col_num, color):
        r = ord(row_char.upper()) - ord('A')
        c = col_num - 1
        if 0 <= r < HEIGHT and 0 <= c < WIDTH:
            board[r][c] = color

    # Setup the initial board state from the diagram
    initial_setup = {
        'B': [('F', 9), ('G', 9), ('H', 9), ('I', 9), ('I', 13), ('J', 8), ('J', 9),
              ('J', 15), ('M', 9), ('N', 9), ('O', 10), ('O', 11)],
        'W': [('I', 10), ('I', 11), ('I', 12), ('J', 10), ('J', 11), ('J', 12), ('J', 13),
              ('J', 14), ('K', 8), ('K', 9), ('K', 10), ('K', 11), ('K', 12), ('L', 8),
              ('L', 9), ('L', 10), ('L', 11), ('L', 12), ('M', 10), ('M', 11), ('M', 12),
              ('N', 10), ('N', 11), ('N', 12), ('O', 12)]
    }

    for color, positions in initial_setup.items():
        for r_char, c_num in positions:
            set_disk(r_char, c_num, color)

    # All 8 directions to check for flips
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    max_flips = 0

    # Iterate over every square on the board
    for r in range(HEIGHT):
        for c in range(WIDTH):
            # We are looking for an empty square to place a black disk
            if board[r][c] == '.':
                current_move_flips = 0
                
                # Check all 8 directions from this empty square
                for dr, dc in directions:
                    flips_in_this_line = []
                    curr_r, curr_c = r + dr, c + dc

                    # Move in the current direction as long as we see opponent's disks
                    while 0 <= curr_r < HEIGHT and 0 <= curr_c < WIDTH and board[curr_r][curr_c] == OPPONENT:
                        flips_in_this_line.append((curr_r, curr_c))
                        curr_r += dr
                        curr_c += dc

                    # If the line ends with one of our own disks, the move is valid in this direction
                    if 0 <= curr_r < HEIGHT and 0 <= curr_c < WIDTH and board[curr_r][curr_c] == PLAYER:
                        current_move_flips += len(flips_in_this_line)
                
                # Update the maximum number of flips found so far
                if current_move_flips > max_flips:
                    max_flips = current_move_flips

    print("The largest number of disks that black can flip in one turn is:")
    print(max_flips)

solve_reversi()