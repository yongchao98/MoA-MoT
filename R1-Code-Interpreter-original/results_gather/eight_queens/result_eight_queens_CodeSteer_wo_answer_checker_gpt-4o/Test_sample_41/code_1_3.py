def solve_n_queens(board):
    n = 8
    queens = []
    restricted = set()
    occupied_cols = set()
    occupied_diag1 = set()  # r - c
    occupied_diag2 = set()  # r + c
    
    # Parse the board to find existing queens and restricted positions
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                queens.append((r, c))
                occupied_cols.add(c)
                occupied_diag1.add(r - c)
                occupied_diag2.add(r + c)
            elif board[r][c] == 'X':
                restricted.add((r, c))
    
    # Backtracking function to place queens
    def place_queens(row):
        if row == n:
            return True
        for col in range(n):
            if (row, col) not in restricted and col not in occupied_cols and (row - col) not in occupied_diag1 and (row + col) not in occupied_diag2:
                # Place the queen
                queens.append((row, col))
                occupied_cols.add(col)
                occupied_diag1.add(row - col)
                occupied_diag2.add(row + col)
                
                # Recurse to place the next queen
                if place_queens(row + 1):
                    return True
                
                # Backtrack
                queens.pop()
                occupied_cols.remove(col)
                occupied_diag1.remove(row - col)
                occupied_diag2.remove(row + col)
        return False
    
    # Start placing queens from the first row
    place_queens(0)
    
    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Given board
board = [
    "XX000000",
    "0000010X",
    "00000000",
    "00000000",
    "00000001",
    "00000000",
    "00000000",
    "00X00000"
]

# Solve the problem
print(solve_n_queens(board))