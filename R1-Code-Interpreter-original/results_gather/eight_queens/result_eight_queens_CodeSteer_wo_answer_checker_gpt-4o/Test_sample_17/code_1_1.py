def solve_n_queens():
    n = 8
    board = [[0] * n for _ in range(n)]
    board[3][3] = 1  # Pre-placed queen

    # Marking the row, column, and diagonals occupied by the pre-placed queen
    rows = [False] * n
    cols = [False] * n
    main_diag = [False] * (2 * n - 1)
    anti_diag = [False] * (2 * n - 1)

    # Mark the pre-placed queen's position
    rows[3] = True
    cols[3] = True
    main_diag[3 - 3 + n - 1] = True
    anti_diag[3 + 3] = True

    # List to store the positions of queens
    positions = [(3, 3)]

    def place_queens(row):
        if row == n:
            return True
        if row == 3:  # Skip the row with the pre-placed queen
            return place_queens(row + 1)
        for col in range(n):
            if not (cols[col] or main_diag[row - col + n - 1] or anti_diag[row + col]):
                # Place the queen
                rows[row] = cols[col] = main_diag[row - col + n - 1] = anti_diag[row + col] = True
                positions.append((row, col))
                print(f"Placing queen at: {row}, {col}")  # Debug statement
                if place_queens(row + 1):
                    return True
                # Backtrack
                rows[row] = cols[col] = main_diag[row - col + n - 1] = anti_diag[row + col] = False
                positions.pop()
                print(f"Backtracking from: {row}, {col}")  # Debug statement
        return False

    # Start placing queens from row 0
    place_queens(0)

    # Output the positions
    print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")

solve_n_queens()