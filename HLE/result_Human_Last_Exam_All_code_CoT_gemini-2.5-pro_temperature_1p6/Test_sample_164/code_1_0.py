def solve_reversi():
    """
    Calculates the largest number of disks black can flip in one turn on the given Reversi board.
    """

    # Helper function to parse coordinates like 'A1', 'T20' into 0-indexed (row, col) tuples.
    def parse_coords(coord_str):
        row_char = coord_str[0]
        col_num = int(coord_str[1:])
        row = ord(row_char) - ord('A')
        col = col_num - 1
        return row, col

    # 1. Digitize the Board
    BOARD_SIZE = 20
    board = [['.' for _ in range(BOARD_SIZE)] for _ in range(BOARD_SIZE)]

    black_pieces_str = "F9 G9 H9 I9 I13 J8 J9 J15 M9 N9 O10 O11".split()
    white_pieces_str = ("I10 I11 I12 J10 J11 J12 J13 J14 K8 K9 K10 K11 K12 L8 L9 L10 "
                        "L11 L12 M10 M11 M12 N10 N11 N12 O12").split()

    for coord in black_pieces_str:
        r, c = parse_coords(coord)
        board[r][c] = 'B'

    for coord in white_pieces_str:
        r, c = parse_coords(coord)
        board[r][c] = 'W'

    # Define the 8 directions to check from any given square.
    directions = [(-1, -1), (-1, 0), (-1, 1),
                  (0, -1),           (0, 1),
                  (1, -1),  (1, 0),  (1, 1)]

    max_flips = 0

    # 2. Identify Potential Moves (iterate over all empty squares)
    for r in range(BOARD_SIZE):
        for c in range(BOARD_SIZE):
            if board[r][c] == '.':
                current_move_flips = 0
                
                # 3. Calculate Flips for Each Move
                for dr, dc in directions:
                    flips_in_this_direction = []
                    
                    line_r, line_c = r + dr, c + dc
                    
                    # Trace a line in the current direction.
                    while 0 <= line_r < BOARD_SIZE and 0 <= line_c < BOARD_SIZE:
                        # If we find a White disk, it's a candidate for flipping.
                        if board[line_r][line_c] == 'W':
                            flips_in_this_direction.append((line_r, line_c))
                        # If we find a Black disk, the sandwich is complete.
                        elif board[line_r][line_c] == 'B':
                            current_move_flips += len(flips_in_this_direction)
                            break 
                        # If we find an empty square, the sandwich is broken.
                        else: # board[line_r][line_c] == '.'
                            break
                        
                        line_r += dr
                        line_c += dc
                
                # 4. Count and Compare
                if current_move_flips > max_flips:
                    max_flips = current_move_flips
    
    # 5. Determine the Maximum
    print(f"The largest number of disks that black can flip is:")
    # The problem asks for the numbers in the final equation, but as there's only one resulting
    # number, we will print it directly.
    print(max_flips)

solve_reversi()