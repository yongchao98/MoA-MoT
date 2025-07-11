import itertools

def count_non_attacking_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5
    toroidal chessboard.
    """
    board_size = 5
    num_queens = 4

    # Create a list of all possible squares (r, c) on the board
    all_squares = [(r, c) for r in range(board_size) for c in range(board_size)]

    # Generate all possible combinations of placing the queens on the board
    possible_placements = itertools.combinations(all_squares, num_queens)

    valid_placement_count = 0
    # Iterate through each possible placement
    for placement in possible_placements:
        # For a placement to be valid, all queens must be on different...
        # 1. rows
        # 2. columns
        # 3. main diagonals (r - c)
        # 4. anti-diagonals (r + c)
        # We use the modulo operator to handle the toroidal "wrap-around" effect.

        rows = {q[0] for q in placement}
        cols = {q[1] for q in placement}
        main_diagonals = {(q[0] - q[1]) % board_size for q in placement}
        anti_diagonals = {(q[0] + q[1]) % board_size for q in placement}

        # A valid placement has no two queens on the same line.
        # This means the number of unique rows, cols, and diagonals
        # must be equal to the number of queens.
        if (len(rows) == num_queens and
            len(cols) == num_queens and
            len(main_diagonals) == num_queens and
            len(anti_diagonals) == num_queens):
            valid_placement_count += 1

    print(f"The number of ways to place {num_queens} non-attacking queens on a {board_size}x{board_size} toroidal chessboard is:")
    print(valid_placement_count)

if __name__ == '__main__':
    count_non_attacking_queens()