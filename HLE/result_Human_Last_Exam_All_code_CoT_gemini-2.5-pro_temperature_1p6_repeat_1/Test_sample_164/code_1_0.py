def solve_reversi():
    """
    Calculates the largest number of disks Black can flip in one turn on the given board.
    """
    # Step 1: Represent the board based on the provided diagram.
    # Rows are A-T (0-19), Columns are 1-20 (0-19).
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
        "........WWWWW.......",  # K
        "........WWWWW.......",  # L
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
    height, width = 20, 20
    player, opponent = 'B', 'W'
    
    # Directions for checking adjacent squares (N, NE, E, SE, S, SW, W, NW)
    directions = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]

    max_flips = 0
    best_move_flips_list = []

    # Step 2: Iterate through every empty square to find potential moves.
    for r in range(height):
        for c in range(width):
            if board[r][c] == '.':
                current_move_total_flips = 0
                current_move_flips_list = []
                
                # Step 3: For each potential move, calculate the number of flips.
                for dr, dc in directions:
                    flips_in_this_direction = []
                    # Move one step in the current direction.
                    line_r, line_c = r + dr, c + dc
                    
                    # Follow the line of opponent's disks.
                    while 0 <= line_r < height and 0 <= line_c < width:
                        if board[line_r][line_c] == opponent:
                            flips_in_this_direction.append((line_r, line_c))
                        elif board[line_r][line_c] == player:
                            # Found an anchor piece. The flips are valid.
                            if flips_in_this_direction:
                                num_flips = len(flips_in_this_direction)
                                current_move_total_flips += num_flips
                                current_move_flips_list.append(num_flips)
                            break
                        else:  # Hit an empty square or invalid piece
                            break
                        line_r += dr
                        line_c += dc
                
                # Step 4: Keep track of the maximum flips found so far.
                if current_move_total_flips > max_flips:
                    max_flips = current_move_total_flips
                    best_move_flips_list = current_move_flips_list

    # Step 5: Format and print the final result.
    if max_flips > 0:
        # Sort for a consistent, readable output equation.
        best_move_flips_list.sort(reverse=True)
        equation_parts = [str(f) for f in best_move_flips_list]
        equation_str = " + ".join(equation_parts)
        print(f"The largest number of disks Black can flip is calculated as:")
        print(f"{equation_str} = {max_flips}")
    else:
        print("0")

solve_reversi()