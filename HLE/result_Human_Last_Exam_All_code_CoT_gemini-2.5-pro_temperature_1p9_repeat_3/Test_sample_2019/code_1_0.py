def solve_toroidal_queens():
    """
    This script calculates the number of ways to place 4 non-attacking queens
    on a 5x5 toroidal chessboard.
    """
    
    # Board dimensions and number of queens
    N = 5
    NUM_QUEENS = 4
    
    # Using a list to hold the count, so it's mutable across recursive calls
    solution_count = [0]

    def is_safe(r, c, queens):
        """
        Checks if placing a queen at (r, c) is safe with respect to the
        list of queens already placed on the toroidal board.
        """
        for qr, qc in queens:
            # Check for same row or column (standard attack)
            if r == qr or c == qc:
                return False
            # Check for toroidal diagonals
            if (r - c) % N == (qr - qc) % N:
                return False
            if (r + c) % N == (qr + qc) % N:
                return False
        return True

    def find_placements(k, start_idx, queens):
        """
        Recursively finds valid placements using backtracking.
        k: number of queens currently placed.
        start_idx: linear board index to start searching from (avoids duplicates).
        queens: a list of (row, col) tuples for placed queens.
        """
        # Base case: if all queens are placed, we found a solution.
        if k == NUM_QUEENS:
            solution_count[0] += 1
            return

        # Iterate through all possible squares from the start index.
        # The board is treated as a 1D array of size N*N.
        for i in range(start_idx, N * N):
            r = i // N
            c = i % N

            # Check if placing a queen at (r, c) is safe.
            if is_safe(r, c, queens):
                # If safe, add the queen to the list of placements.
                queens.append((r, c))
                # Recurse to place the next queen.
                # Start the search for the next queen from the next index (i + 1).
                find_placements(k + 1, i + 1, queens)
                # Backtrack: remove the queen to explore other configurations.
                queens.pop()

    # Initial call to start the backtracking process.
    # Start with 0 queens placed, from index 0, with an empty list of queens.
    find_placements(0, 0, [])

    # The problem asks for the total number of ways.
    # The variable solution_count[0] holds this number.
    print(solution_count[0])

solve_toroidal_queens()
<<<200>>>