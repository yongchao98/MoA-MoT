import collections

def solve_reversi():
    """
    Solves the Reversi problem by finding the best move for black.
    """
    board_size = 20
    board = [['.' for _ in range(board_size)] for _ in range(board_size)]
    player = 'B'
    opponent = 'W'

    # Board setup from the problem description
    # Mapping 'A'-'T' to 0-19 and '1'-'20' to 0-19
    pieces = {
        'F9': 'B', 'G9': 'B', 'H9': 'B', 'I9': 'B', 'I10': 'W', 'I11': 'W',
        'I12': 'W', 'I13': 'B', 'J8': 'B', 'J9': 'B', 'J10': 'W', 'J11': 'W',
        'J12': 'W', 'J13': 'W', 'J14': 'W', 'J15': 'B', 'K8': 'W', 'K9': 'W',
        'K10': 'W', 'K11': 'W', 'K12': 'W', 'L8': 'W', 'L9': 'W', 'L10': 'W',
        'L11': 'W', 'L12': 'W', 'M9': 'B', 'M10': 'W', 'M11': 'W', 'M12': 'W',
        'N9': 'B', 'N10': 'W', 'N11': 'W', 'N12': 'W', 'O10': 'B', 'O11': 'B', 'O12': 'W'
    }

    for pos, piece in pieces.items():
        row_char = pos[0]
        col_str = pos[1:]
        r = ord(row_char) - ord('A')
        c = int(col_str) - 1
        board[r][c] = piece

    # Directions: (dr, dc)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1),
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]

    max_flips = 0
    best_move_details = []

    # Iterate through all empty cells
    for r in range(board_size):
        for c in range(board_size):
            if board[r][c] == '.':
                current_flips = 0
                move_details = []
                
                # Check all 8 directions for this potential move
                for dr, dc in directions:
                    flips_in_direction = 0
                    line = []
                    
                    # Trace the line
                    nr, nc = r + dr, c + dc
                    while 0 <= nr < board_size and 0 <= nc < board_size:
                        if board[nr][nc] == opponent:
                            line.append((nr, nc))
                        elif board[nr][nc] == player:
                            flips_in_direction = len(line)
                            break
                        else: # Empty or off-board
                            break
                        nr, nc = nr + dr, nc + dc

                    if flips_in_direction > 0:
                        current_flips += flips_in_direction
                        move_details.append(flips_in_direction)
                
                # Check if this move is the best so far
                if current_flips > max_flips:
                    max_flips = current_flips
                    best_move_details = sorted(move_details, reverse=True)

    # Output the results
    if max_flips > 0:
        equation_str = " + ".join(map(str, best_move_details))
        print(f"The largest number of disks black can flip is {max_flips}.")
        print(f"This is achieved by flipping disks in several directions, resulting in the sum:")
        print(f"{equation_str} = {max_flips}")
    else:
        print("There are no valid moves for black.")

solve_reversi()