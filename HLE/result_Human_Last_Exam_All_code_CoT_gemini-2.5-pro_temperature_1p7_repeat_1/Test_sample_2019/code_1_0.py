def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4
    
    # Using a list to hold the count so it can be modified within the nested function
    solution_count = [0]

    def is_attacked(p1, p2):
        """
        Checks if two queens at positions p1 and p2 attack each other on the NxN toroidal board.
        p1 and p2 are tuples (row, col).
        """
        r1, c1 = p1
        r2, c2 = p2

        # Check for row or column attack
        if r1 == r2 or c1 == c2:
            return True

        # Check for toroidal diagonal attack
        dr = abs(r1 - r2)
        dc = abs(c1 - c2)
        if min(dr, N - dr) == min(dc, N - dc):
            return True

        return False

    def backtrack(k, start_index, placements):
        """
        Recursively finds valid placements for the queens.
        k: The number of queens still to be placed.
        start_index: The linear index on the board (0 to 24) to start searching from.
                     This prevents duplicate combinations.
        placements: A list of (row, col) tuples of queens already placed.
        """
        # Base case: If all queens have been placed, we found a solution.
        if k == 0:
            solution_count[0] += 1
            return

        # Optimization: if not enough squares are left for the remaining queens, prune this branch.
        if N * N - start_index < k:
            return

        # Iterate through all possible cells starting from start_index
        for i in range(start_index, N * N):
            r = i // N
            c = i % N
            current_pos = (r, c)
            
            # Check if placing a queen at current_pos is safe
            is_safe = True
            for queen_pos in placements:
                if is_attacked(current_pos, queen_pos):
                    is_safe = False
                    break
            
            # If it's safe, place the queen and recurse
            if is_safe:
                placements.append(current_pos)
                backtrack(k - 1, i + 1, placements)
                # Backtrack: remove the queen to explore other possibilities
                placements.pop()

    # Start the backtracking process to place NUM_QUEENS queens
    backtrack(NUM_QUEENS, 0, [])

    print(f"Number of ways to place {NUM_QUEENS} non-attacking queens on a {N}x{N} toroidal chessboard: {solution_count[0]}")

solve_toroidal_queens()
<<<400>>>