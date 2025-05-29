def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row, n), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
            
    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row, n), range(col, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens():
    n = 8
    # Initialize board with given constraints
    board = [
        [0,0,0,-1,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,1,0,-1,-1,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    def solve_util(col):
        if col >= n:
            return True
            
        for row in range(n):
            if board[row][col] != -1 and board[row][col] != 1:
                # Skip if there's already a queen in this row
                row_has_queen = False
                for j in range(n):
                    if board[row][j] == 1:
                        row_has_queen = True
                        break
                if row_has_queen:
                    continue
                
                if is_safe(board, row, col, n):
                    board[row][col] = 1
                    if solve_util(col + 1):
                        return True
                    board[row][col] = 0
        return False
    
    # Find first empty column
    start_col = 0
    queens_placed = []
    
    # First, record the pre-placed queen
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens_placed.append(f"{i} {j}")
    
    if solve_util(start_col):
        # Get positions of queens
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1 and f"{i} {j}" not in queens_placed:
                    queens_placed.append(f"{i} {j}")
        queens_placed.sort(key=lambda x: (int(x.split()[0]), int(x.split()[1])))
        print(",".join(queens_placed))
    else:
        print("No solution exists")

solve_n_queens()