import itertools

def count_toroidal_queens_placements():
    """
    Calculates the number of ways to place a given number of non-attacking queens
    on a toroidal chessboard of a given size.
    """
    board_size = 5
    num_queens = 4

    # Create a list of all possible square coordinates, e.g., (0, 0), (0, 1), ...
    all_squares = [(r, c) for r in range(board_size) for c in range(board_size)]

    valid_placements_count = 0

    # Generate all unique combinations of squares to place the queens
    # This is equivalent to "25 choose 4" for a 5x5 board and 4 queens.
    for queen_positions in itertools.combinations(all_squares, num_queens):
        is_valid = True
        
        # Check for attacks between every pair of queens in the current placement
        # For 4 queens, there are C(4, 2) = 6 pairs to check.
        for q1, q2 in itertools.combinations(queen_positions, 2):
            r1, c1 = q1
            r2, c2 = q2

            # Check for row or column attacks
            if r1 == r2 or c1 == c2:
                is_valid = False
                break
            
            # Check for toroidal diagonal attacks using modular arithmetic
            if (r1 - c1) % board_size == (r2 - c2) % board_size:
                is_valid = False
                break
            
            if (r1 + c1) % board_size == (r2 + c2) % board_size:
                is_valid = False
                break

        if is_valid:
            valid_placements_count += 1

    # Print the final result in a descriptive sentence
    print(f"For a {board_size}x{board_size} toroidal board, the number of ways to place {num_queens} non-attacking queens is {valid_placements_count}.")

if __name__ == '__main__':
    count_toroidal_queens_placements()