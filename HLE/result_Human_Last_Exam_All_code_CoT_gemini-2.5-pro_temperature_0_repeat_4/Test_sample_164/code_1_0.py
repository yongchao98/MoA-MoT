def solve_reversi():
    """
    Calculates the largest number of disks black can flip in one turn on a custom Reversi board.
    """
    # 1. Represent the Board
    # Board size is 20x20, rows A-T (0-19), cols 1-20 (0-19)
    BOARD_SIZE = 20
    board = [['.' for _ in range(BOARD_SIZE)] for _ in range(BOARD_SIZE)]

    # The given positions of the pieces
    pieces = {
        ('F', 9): 'B', ('G', 9): 'B', ('H', 9): 'B', ('I', 9): 'B',
        ('I', 10): 'W', ('I', 11): 'W', ('I', 12): 'W', ('I', 13): 'B',
        ('J', 8): 'B', ('J', 9): 'B', ('J', 10): 'W', ('J', 11): 'W',
        ('J', 12): 'W', ('J', 13): 'W', ('J', 14): 'W', ('J', 15): 'B',
        ('K', 8): 'W', ('K', 9): 'W', ('K', 10): 'W', ('K', 11): 'W', ('K', 12): 'W',
        ('L', 8): 'W', ('L', 9): 'W', ('L', 10): 'W', ('L', 11): 'W', ('L', 12): 'W',
        ('M', 9): 'B', ('M', 10): 'W', ('M', 11): 'W', ('M', 12): 'W',
        ('N', 9): 'B', ('N', 10): 'W', ('N', 11): 'W', ('N', 12): 'W',
        ('O', 10): 'B', ('O', 11): 'B', ('O', 12): 'W'
    }

    # Populate the board grid
    for (row_char, col_num), color in pieces.items():
        row_idx = ord(row_char) - ord('A')
        col_idx = col_num - 1
        if 0 <= row_idx < BOARD_SIZE and 0 <= col_idx < BOARD_SIZE:
            board[row_idx][col_idx] = color

    max_flips = 0
    best_breakdown = []
    
    # Define the 8 directions to check
    directions = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
    
    # 2. Identify All Possible Moves & 3. Calculate Flips
    for r in range(BOARD_SIZE):
        for c in range(BOARD_SIZE):
            # Consider placing a black disk on an empty square
            if board[r][c] == '.':
                current_flips_total = 0
                current_breakdown = []
                
                for dr, dc in directions:
                    flips_in_direction = 0
                    line_to_check = []
                    
                    curr_r, curr_c = r + dr, c + dc
                    
                    # Move along the current direction
                    while 0 <= curr_r < BOARD_SIZE and 0 <= curr_c < BOARD_SIZE:
                        if board[curr_r][curr_c] == 'W':
                            line_to_check.append((curr_r, curr_c))
                        elif board[curr_r][curr_c] == 'B':
                            # Found an anchor, this is a valid flip
                            if line_to_check:
                                flips_in_direction = len(line_to_check)
                                current_flips_total += flips_in_direction
                                current_breakdown.append(flips_in_direction)
                            break # Stop searching in this direction
                        else: # Found an empty square or went off board
                            break # Stop searching in this direction
                            
                        curr_r, curr_c = curr_r + dr, curr_c + dc
                
                # 4. Find the Maximum
                if current_flips_total > max_flips:
                    max_flips = current_flips_total
                    best_breakdown = current_breakdown

    # 5. Output the Result
    if max_flips > 0:
        # Sort for consistent output, largest numbers first
        best_breakdown.sort(reverse=True)
        equation_parts = [str(x) for x in best_breakdown]
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {max_flips}")
    else:
        print("There are no valid moves for black to flip any disks.")

solve_reversi()