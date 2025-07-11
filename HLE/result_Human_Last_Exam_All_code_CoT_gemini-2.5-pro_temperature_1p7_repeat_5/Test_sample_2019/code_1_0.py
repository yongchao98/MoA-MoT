def solve_toroidal_queens():
    """
    Calculates the number of ways to place K non-attacking queens on an N x N toroidal board.
    """
    N = 5
    K = 4
    
    # Using a list to hold the count so it can be modified within the nested function.
    solution_count = [0]

    def is_safe(r_new, c_new, placements):
        """
        Checks if placing a queen at (r_new, c_new) is safe with respect to existing placements.
        """
        for r_old, c_old in placements:
            # Check for row or column attack
            if r_new == r_old or c_new == c_old:
                return False

            # Check for toroidal diagonal attack
            d_row = abs(r_new - r_old)
            d_col = abs(c_new - c_old)
            
            toroidal_d_row = min(d_row, N - d_row)
            toroidal_d_col = min(d_col, N - d_col)

            if toroidal_d_row == toroidal_d_col:
                return False
        return True

    def find_placements(start_index, placements):
        """
        A recursive backtracking function to find all valid placements.
        - start_index: The linear index on the board (0-24) from which to start searching.
        - placements: A list of (row, col) tuples for queens already on the board.
        """
        if len(placements) == K:
            solution_count[0] += 1
            return

        # Iterate through the remaining squares on the board
        for i in range(start_index, N * N):
            r = i // N
            c = i % N

            if is_safe(r, c, placements):
                # If safe, place the queen
                placements.append((r, c))
                # Recurse to find placements for the next queen
                find_placements(i + 1, placements)
                # Backtrack: remove the queen to explore other possibilities
                placements.pop()

    # Start the search from the first square (index 0) with an empty set of placements.
    find_placements(0, [])
    
    final_count = solution_count[0]
    print(f"On a {N}x{N} toroidal chessboard, there are {final_count} ways to place {K} non-attacking queens.")

solve_toroidal_queens()
<<<200>>>