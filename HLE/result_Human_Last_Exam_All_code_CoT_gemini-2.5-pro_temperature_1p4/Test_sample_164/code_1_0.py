def solve_reversi():
    """
    Analyzes a Reversi board to find the best move for the black player.
    """
    # 1. Represent the Board
    rows, cols = 20, 20
    board = [['.' for _ in range(cols)] for _ in range(rows)]
    player = 'B'
    opponent = 'W'

    # Helper to convert board coordinates like 'A', 'B', etc., to row indices 0, 1, etc.
    def char_to_row(c):
        return ord(c) - ord('A')

    # The initial positions of the disks from the diagram.
    # Coordinates are (row_char, col_num, player_color).
    piece_data = [
        ('F', 9, 'B'), ('G', 9, 'B'), ('H', 9, 'B'), ('I', 9, 'B'),
        ('I', 10, 'W'), ('I', 11, 'W'), ('I', 12, 'W'), ('I', 13, 'B'),
        ('J', 8, 'B'), ('J', 9, 'B'), ('J', 10, 'W'), ('J', 11, 'W'),
        ('J', 12, 'W'), ('J', 13, 'W'), ('J', 14, 'W'), ('J', 15, 'B'),
        ('K', 8, 'W'), ('K', 9, 'W'), ('K', 10, 'W'), ('K', 11, 'W'),
        ('K', 12, 'W'),
        ('L', 8, 'W'), ('L', 9, 'W'), ('L', 10, 'W'), ('L', 11, 'W'),
        ('L', 12, 'W'),
        ('M', 9, 'B'), ('M', 10, 'W'), ('M', 11, 'W'), ('M', 12, 'W'),
        ('N', 9, 'B'), ('N', 10, 'W'), ('N', 11, 'W'), ('N', 12, 'W'),
        ('O', 10, 'B'), ('O', 11, 'B'), ('O', 12, 'W')
    ]

    for r_char, c_num, color in piece_data:
        r, c = char_to_row(r_char), c_num - 1
        if 0 <= r < rows and 0 <= c < cols:
            board[r][c] = color

    # Directions to check: N, S, E, W, NE, NW, SE, SW
    directions = [(-1, 0), (1, 0), (0, 1), (0, -1),
                  (-1, 1), (-1, -1), (1, 1), (1, -1)]

    max_flips = 0
    best_move_details = []

    # 2. Identify All Possible Moves and 3. Calculate Flips
    for r in range(rows):
        for c in range(cols):
            # A move is only possible on an empty square
            if board[r][c] == '.':
                flips_from_this_move = []
                
                # Check all 8 directions from the current empty square
                for dr, dc in directions:
                    discs_to_flip = []
                    r_check, c_check = r + dr, c + dc

                    # Find a continuous line of opponent's disks
                    while 0 <= r_check < rows and 0 <= c_check < cols and board[r_check][c_check] == opponent:
                        discs_to_flip.append((r_check, c_check))
                        r_check += dr
                        c_check += dc
                    
                    # If the line is anchored by one of the player's disks, it's a valid flip
                    if (0 <= r_check < rows and 0 <= c_check < cols and 
                        board[r_check][c_check] == player and discs_to_flip):
                        flips_from_this_move.append(len(discs_to_flip))
                
                # 4. Find the Maximum
                current_total_flips = sum(flips_from_this_move)
                if current_total_flips > max_flips:
                    max_flips = current_total_flips
                    best_move_details = flips_from_this_move

    # 5. Output the Result
    # Sort the flips per direction for a consistent output format
    equation_numbers = sorted(best_move_details, reverse=True)
    
    print("The final equation is: ", end="")
    for i, num in enumerate(equation_numbers):
        print(num, end="")
        if i < len(equation_numbers) - 1:
            print(" + ", end="")
    
    print(f" = {max_flips}")
    print(f"\nThe largest number of disks that can be flipped is {max_flips}.")


solve_reversi()