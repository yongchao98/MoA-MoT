def is_safe(board, row, col, n, cols, diag1, diag2):
    # Check if position is marked with 'X'
    if board[row][col] == 'X':
        return False
    
    # Check if the column is occupied
    if col in cols:
        return False
    
    # Check diagonals
    if (row + col) in diag1 or (row - col) in diag2:
        return False
    
    return True

def solve_queens(board, row, n, cols, diag1, diag2, queens):
    if len(queens) == 8:
        return True
    
    if row >= n:
        return False
    
    # If this row already has a queen, move to next row
    for i in range(n):
        if board[row][i] == '1':
            queens.add((row, i))
            cols.add(i)
            diag1.add(row + i)
            diag2.add(row - i)
            return solve_queens(board, row + 1, n, cols, diag1, diag2, queens)
    
    # Try placing queen in each column of this row
    for col in range(n):
        if is_safe(board, row, col, n, cols, diag1, diag2):
            # Place queen
            board[row][col] = '1'
            queens.add((row, col))
            cols.add(col)
            diag1.add(row + col)
            diag2.add(row - col)
            
            if solve_queens(board, row + 1, n, cols, diag1, diag2, queens):
                return True
            
            # Remove queen if solution not found
            board[row][col] = '0'
            queens.remove((row, col))
            cols.remove(col)
            diag1.remove(row + col)
            diag2.remove(row - col)
    
    return False

# Initialize board
board = [
    list('00X00000'),
    list('01000000'),
    list('00000000'),
    list('0000000X'),
    list('00000000'),
    list('00000000'),
    list('00001000'),
    list('00000000')
]

n = 8
cols = set()
diag1 = set()  # row + col
diag2 = set()  # row - col
queens = set()

# Add pre-placed queens
if board[1][1] == '1':
    queens.add((1, 1))
    cols.add(1)
    diag1.add(2)
    diag2.add(0)

if board[6][4] == '1':
    queens.add((6, 4))
    cols.add(4)
    diag1.add(10)
    diag2.add(2)

solve_queens(board, 0, n, cols, diag1, diag2, queens)

# Sort queens by row for consistent output
queen_positions = sorted(list(queens))
result = ', '.join(f"{row} {col}" for row, col in queen_positions)
print(f"<<<{result}>>>")