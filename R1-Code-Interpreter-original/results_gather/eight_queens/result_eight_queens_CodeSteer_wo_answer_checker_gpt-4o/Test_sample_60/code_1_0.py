def solve_n_queens():
    # Initialize the board with the given configuration
    board = [[0 for _ in range(8)] for _ in range(8)]
    board[7][7] = 1  # The last cell is already occupied by a queen

    # Sets to track columns and diagonals that are occupied
    columns = {7}
    pos_diagonals = {7 + 7}  # r + c
    neg_diagonals = {7 - 7}  # r - c

    # List to store the positions of queens
    queen_positions = [(7, 7)]

    def place_queens(row):
        if row == 8:
            return True

        # Skip the row if it already has a queen
        if row == 7:
            return place_queens(row + 1)

        for col in range(8):
            if col in columns or (row + col) in pos_diagonals or (row - col) in neg_diagonals:
                continue

            # Place the queen
            board[row][col] = 1
            columns.add(col)
            pos_diagonals.add(row + col)
            neg_diagonals.add(row - col)
            queen_positions.append((row, col))

            # Debug statement to trace the placement
            print(f"Placing queen at: ({row}, {col})")

            # Recurse to place the next queen
            if place_queens(row + 1):
                return True

            # Backtrack
            board[row][col] = 0
            columns.remove(col)
            pos_diagonals.remove(row + col)
            neg_diagonals.remove(row - col)
            queen_positions.pop()

            # Debug statement to trace backtracking
            print(f"Backtracking from: ({row}, {col})")

        return False

    # Start placing queens from the first row
    place_queens(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queen_positions))
    print(f'<<<{result}>>>')

solve_n_queens()