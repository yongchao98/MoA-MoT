import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4
    solution_count = 0

    # Create a list of all possible cells on the board, e.g., (0, 0), (0, 1), ...
    all_cells = [(r, c) for r in range(N) for c in range(N)]

    # Generate all combinations of placing NUM_QUEENS on the available cells
    # A placement is a tuple of queen positions, e.g., ((r1, c1), (r2, c2), ...)
    for placement in itertools.combinations(all_cells, NUM_QUEENS):
        is_valid = True
        
        # Check every pair of queens in the current placement for attacks
        for queen1, queen2 in itertools.combinations(placement, 2):
            r1, c1 = queen1
            r2, c2 = queen2

            # Check for row, column, and toroidal diagonal/anti-diagonal attacks
            if (r1 == r2 or
                    c1 == c2 or
                    (r1 - c1) % N == (r2 - c2) % N or
                    (r1 + c1) % N == (r2 + c2) % N):
                is_valid = False
                break  # An attack was found, no need to check other pairs in this placement

        if is_valid:
            solution_count += 1
            
    # The final equation shows the parameters of the problem and its solution
    print(f"Ways to place {NUM_QUEENS} queens on a {N}x{N} toroidal board = {solution_count}")

solve_toroidal_queens()
