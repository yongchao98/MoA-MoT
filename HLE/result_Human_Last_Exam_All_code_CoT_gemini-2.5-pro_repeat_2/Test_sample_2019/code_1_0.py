def solve_n_queens_toroidal():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    board_size = 5
    num_queens = 4

    def is_safe(queen_pos, placed_queens):
        """
        Checks if a queen at queen_pos is safe from the queens in placed_queens.
        """
        r2, c2 = queen_pos
        for r1, c1 in placed_queens:
            # Check for row or column attack
            if r1 == r2 or c1 == c2:
                return False
            # Check for toroidal diagonal attacks
            if (r1 - c1) % board_size == (r2 - c2) % board_size:
                return False
            if (r1 + c1) % board_size == (r2 + c2) % board_size:
                return False
        return True

    def count_placements(queens_to_place, start_index, placed_queens):
        """
        Recursively counts valid queen placements using backtracking.
        """
        # Base case: if all queens have been placed, we found one valid solution.
        if queens_to_place == 0:
            return 1

        count = 0
        # Iterate over possible squares, starting from start_index to avoid duplicates.
        for i in range(start_index, board_size * board_size):
            # Optimization: if remaining squares are not enough for remaining queens, stop.
            if (board_size * board_size - i) < queens_to_place:
                break
                
            r, c = i // board_size, i % board_size

            if is_safe((r, c), placed_queens):
                # Place the queen
                placed_queens.append((r, c))
                # Recurse to place the next queen
                count += count_placements(queens_to_place - 1, i + 1, placed_queens)
                # Backtrack: remove the queen to explore other options
                placed_queens.pop()
        
        return count

    total_ways = count_placements(num_queens, 0, [])
    
    # The final output explains the result, including the numbers from the problem statement.
    print(f"On a {board_size}x{board_size} toroidal chessboard, there are {total_ways} ways to place {num_queens} non-attacking queens.")

if __name__ == '__main__':
    solve_n_queens_toroidal()