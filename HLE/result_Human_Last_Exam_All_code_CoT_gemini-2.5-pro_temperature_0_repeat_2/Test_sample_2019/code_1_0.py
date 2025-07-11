import itertools

def count_non_attacking_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5
    toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4

    # Create a list of all possible squares on the board, e.g., (0, 0), (0, 1), ...
    squares = [(r, c) for r in range(N) for c in range(N)]

    # Generate all unique combinations of placing NUM_QUEENS on the squares
    possible_placements = itertools.combinations(squares, NUM_QUEENS)

    valid_placement_count = 0

    # Iterate through each possible placement
    for placement in possible_placements:
        is_valid = True
        # Check for attacks between any pair of queens in the current placement.
        # We can use itertools.combinations again to get all unique pairs of queens.
        for q1, q2 in itertools.combinations(placement, 2):
            r1, c1 = q1
            r2, c2 = q2

            # Check for attacks: row, column, and the two toroidal diagonals
            if (r1 == r2 or
                    c1 == c2 or
                    (r1 - c1) % N == (r2 - c2) % N or
                    (r1 + c1) % N == (r2 + c2) % N):
                is_valid = False
                break  # An attack was found, so this placement is invalid.

        if is_valid:
            valid_placement_count += 1

    # Print the final result including the parameters of the problem
    print(f"On a {N}x{N} toroidal board, the number of ways to place {NUM_QUEENS} non-attacking queens is:")
    print(valid_placement_count)

if __name__ == "__main__":
    count_non_attacking_queens()