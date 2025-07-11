import collections

def solve_reversi_max_flips():
    """
    This script calculates the maximum number of disks the black player can flip in a single turn
    on the given Reversi board configuration.

    It works by:
    1. Setting up a 2D list to represent the 20x20 game board.
    2. Populating the board with Black ('B') and White ('W') disks according to the problem description.
    3. Iterating through every empty ('.') square on the board to evaluate it as a potential move for Black.
    4. For each potential move, checking all 8 directions to see if any lines of White disks are
       outflanked and counting the number of disks that would be flipped.
    5. Keeping track of the move that yields the highest number of flips.
    6. Finally, printing the maximum number of flips found, the coordinates of the move,
       and the equation that shows how the total is calculated from flips in each direction.
    """
    # Board representation: 20x20 grid, rows 0-19 (A-T), cols 0-19 (1-20)
    board = [['.' for _ in range(20)] for _ in range(20)]
    player = 'B'
    opponent = 'W'

    # Place Black disks ('B') based on the diagram
    board[5][8] = 'B'   # F9
    board[6][8] = 'B'   # G9
    board[7][8] = 'B'   # H9
    board[8][8] = 'B'   # I9
    board[8][12] = 'B'  # I13
    board[9][7] = 'B'   # J8
    board[9][8] = 'B'   # J9
    board[9][14] = 'B'  # J15
    board[12][8] = 'B'  # M9
    board[13][8] = 'B'  # N9
    board[14][9] = 'B'  # O10
    board[14][10] = 'B' # O11

    # Place White disks ('W') based on the diagram
    board[8][9] = 'W'   # I10
    board[8][10] = 'W'  # I11
    board[8][11] = 'W'  # I12
    board[9][9] = 'W'   # J10
    board[9][10] = 'W'  # J11
    board[9][11] = 'W'  # J12
    board[9][12] = 'W'  # J13
    board[9][13] = 'W'  # J14
    board[10][8] = 'W'  # K9
    board[10][9] = 'W'  # K10
    board[10][10] = 'W' # K11
    board[10][11] = 'W' # K12
    board[11][8] = 'W'  # L9
    board[11][9] = 'W'  # L10
    board[11][10] = 'W' # L11
    board[11][11] = 'W' # L12
    board[12][9] = 'W'  # M10
    board[12][10] = 'W' # M11
    board[12][11] = 'W' # M12
    board[13][9] = 'W'  # N10
    board[13][10] = 'W' # N11
    board[13][11] = 'W' # N12
    board[14][11] = 'W' # O12

    # Directions represented as (row_change, col_change)
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

    max_flips = 0
    best_move_coords = None
    best_move_breakdown = {}

    # Iterate over every cell to find the best possible move
    for r in range(20):
        for c in range(20):
            if board[r][c] == '.':
                current_total_flips = 0
                current_breakdown = collections.defaultdict(int)

                for dr, dc in directions:
                    line_to_flip = []
                    # Trace a line from the potential move in the current direction
                    curr_r, curr_c = r + dr, c + dc
                    while 0 <= curr_r < 20 and 0 <= curr_c < 20:
                        if board[curr_r][curr_c] == opponent:
                            line_to_flip.append((curr_r, curr_c))
                        elif board[curr_r][curr_c] == player:
                            # A valid flanking line is found
                            if line_to_flip:
                                flips_in_direction = len(line_to_flip)
                                current_breakdown[(dr, dc)] = flips_in_direction
                                current_total_flips += flips_in_direction
                            break
                        else: # Empty square or off-board breaks the line
                            break
                        curr_r += dr
                        curr_c += dc

                if current_total_flips > max_flips:
                    max_flips = current_total_flips
                    best_move_coords = (chr(ord('A') + r), c + 1)
                    best_move_breakdown = current_breakdown

    # Print the final result and the calculation breakdown
    if max_flips > 0:
        row_label, col_label = best_move_coords
        print(f"The largest number of disks that black can flip in one turn is {max_flips}.")
        print(f"This is achieved by placing a black disk at {row_label}{col_label}.")

        breakdown_parts = [str(num) for num in best_move_breakdown.values()]
        equation = " + ".join(breakdown_parts)
        print(f"The calculation for the total flips is: {equation} = {max_flips}")
    else:
        print("There are no legal moves for black to flip any disks.")

solve_reversi_max_flips()
>>>7