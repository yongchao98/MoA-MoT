def solve_4_queens_toroidal():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal board.
    """
    board_size = 5
    num_queens = 4
    
    # Using a list to store the positions (r, c) of placed queens.
    placements = []
    
    # We use a list to hold the final count, as nonlocal is not available in Python 2.
    # This makes the code more portable and is a common pattern for modifying outer-scope integers.
    solution_count = [0]

    def is_safe(new_pos, existing_placements):
        """
        Checks if a new queen at new_pos attacks any existing queens.
        """
        r2, c2 = new_pos
        for r1, c1 in existing_placements:
            # Check row and column conflict
            if r1 == r2 or c1 == c2:
                return False
            # Check toroidal diagonal conflict
            if (r1 - c1) % board_size == (r2 - c2) % board_size:
                return False
            if (r1 + c1) % board_size == (r2 + c2) % board_size:
                return False
        return True

    def backtrack(k, start_square_index):
        """
        Recursively places k queens on the board.
        
        k: Number of queens left to place.
        start_square_index: The board square index (0-24) to start searching from.
                            This prevents duplicate combinations.
        """
        # Base case: If all queens have been placed, we found a solution.
        if k == 0:
            solution_count[0] += 1
            return

        # Optimization: If the number of remaining squares is less than the number of
        # queens we still need to place, we can prune this search branch.
        squares_available = board_size * board_size - start_square_index
        if squares_available < k:
            return

        # Iterate through all possible squares for the current queen.
        for i in range(start_square_index, board_size * board_size):
            r, c = i // board_size, i % board_size
            
            # Check if this position is safe with respect to the queens already placed.
            if is_safe((r, c), placements):
                # Place the queen
                placements.append((r, c))
                # Recursively try to place the next queen, starting from the next square
                backtrack(k - 1, i + 1)
                # Backtrack: remove the queen to explore other possibilities
                placements.pop()

    # Start the backtracking process to place num_queens, starting from square 0.
    backtrack(num_queens, 0)
    
    # The final count is the first element of our list
    final_count = solution_count[0]
    
    # Since the problem doesn't have a simple formula, the result is the computed count.
    # The instruction "output each number in the final equation" isn't directly
    # applicable as we use a search algorithm. So we print the final result.
    print(f"The total number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard is:")
    print(final_count)

solve_4_queens_toroidal()