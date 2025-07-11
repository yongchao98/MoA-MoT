import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 
    toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4
    valid_placements = 0

    # All possible squares on the board, represented by a number from 0 to 24.
    all_squares = range(N * N)

    # Generate all combinations of choosing NUM_QUEENS squares from all_squares.
    # For a 5x5 board and 4 queens, this is C(25, 4).
    queen_combinations = itertools.combinations(all_squares, NUM_QUEENS)

    def is_valid(combination):
        """
        Checks if a given combination of queen positions is valid (non-attacking).
        'combination' is a tuple of linear indices, e.g., (0, 7, 13, 19).
        """
        # Convert linear indices to (row, col) coordinates.
        positions = []
        for pos in combination:
            r = pos // N
            c = pos % N
            positions.append((r, c))

        # Check for attacks between every pair of queens in this combination.
        for i in range(NUM_QUEENS):
            for j in range(i + 1, NUM_QUEENS):
                r1, c1 = positions[i]
                r2, c2 = positions[j]

                # Check for row or column attack.
                if r1 == r2 or c1 == c2:
                    return False

                # Check for toroidal diagonal attack.
                if (r1 - c1) % N == (r2 - c2) % N:
                    return False
                if (r1 + c1) % N == (r2 + c2) % N:
                    return False
        
        # If no attacks are found between any pair, the placement is valid.
        return True

    # Iterate through all combinations and count the valid ones.
    for combo in queen_combinations:
        if is_valid(combo):
            valid_placements += 1

    # The problem asks to output the numbers in the final equation.
    # Here, the calculation is a count, so we print the result in a sentence
    # that includes the key numbers from the problem.
    print(f"On a {N}x{N} toroidal board, there are {valid_placements} ways to place {NUM_QUEENS} non-attacking queens.")

solve_toroidal_queens()