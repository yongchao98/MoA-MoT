def is_safe(board, row, col):
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_queens(board, queens_placed=0):
    if queens_placed == 8:
        positions = []
        for i in range(8):
            for j in range(8):
                if board[i][j] == 1:
                    positions.append((i, j))
        positions.sort()
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
        return True

    # Skip rows that already have queens
    for row in range(8):
        has_queen = False
        for col in range(8):
            if board[row][col] == 1:
                has_queen = True
                break
        if not has_queen:
            for col in range(8):
                if board[row][col] != 'X' and board[row][col] != 1:
                    if is_safe(board, row, col):
                        board[row][col] = 1
                        if solve_queens(board, queens_placed + 1):
                            return True
                        board[row][col] = 0
    return False

# Initialize board with the given queen and blocked position
board = [
    [0,0,0,0,0,1,0,0],  # Queen at (0,5)
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,'X',0,0,0,0,0,0],  # Blocked at (3,1)
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

# Count initial queens
initial_queens = sum(row.count(1) for row in board)
solve_queens(board, initial_queens)