def solve_n_queens_toroidal():
    """
    Calculates the number of ways to place K non-attacking queens on an
    N x N toroidal chessboard using a backtracking algorithm.
    """
    N = 5  # The size of the board (N x N)
    K = 4  # The number of queens to place

    def is_safe(r, c, placements):
        """
        Checks if placing a queen at (r, c) is safe with respect to existing
        placements on the toroidal board.
        """
        for qr, qc in placements:
            # Check for row or column attacks
            if r == qr or c == qc:
                return False

            # Check for toroidal diagonal attacks.
            # (a % N + N) % N is a robust way to compute a modulo n 
            # for both positive and negative a.
            if ((r - c) % N + N) % N == ((qr - qc) % N + N) % N:
                return False
            if (r + c) % N == (qr + qc) % N:
                return False
        return True

    def count_placements_recursive(queens_to_place, start_index, placements):
        """
        Recursively finds the number of valid ways to place the remaining queens.

        :param queens_to_place: The number of queens yet to be placed.
        :param start_index: The board square index (0-24) to start searching from.
                            This avoids counting permutations of the same placement.
        :param placements: A list of (row, col) tuples for queens already on the board.
        :return: The number of valid placements found from the current state.
        """
        # Base case: If all K queens have been placed, we found one valid solution.
        if queens_to_place == 0:
            return 1

        count = 0
        # Iterate through squares, starting from 'start_index' to find a valid
        # spot for the next queen.
        for i in range(start_index, N * N):
            r = i // N
            c = i % N

            # If it's safe to place a queen here...
            if is_safe(r, c, placements):
                # ...place it and recurse to find placements for the rest of the queens.
                new_placements = placements + [(r, c)]
                count += count_placements_recursive(queens_to_place - 1, i + 1, new_placements)
        
        return count

    # Start the calculation
    total_ways = count_placements_recursive(K, 0, [])

    # The problem asks to "output each number in the final equation!". 
    # Since the answer is a single number, this is the final result of our calculation.
    print(total_ways)


solve_n_queens_toroidal()