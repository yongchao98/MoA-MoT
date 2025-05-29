def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with given queens and blocked positions
    board = [
        [0,0,0,0,0,0,0,-1],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [1,0,0,0,0,0,0,0],
        [0,0,0,0,-1,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    queens_positions = []
    # Add existing queens
    queens_positions.append((4,0))
    queens_positions.append((3,4))
    
    def solve_util(board, queens_left):
        if queens_left == 0:
            return True
            
        for i in range(n):
            for j in range(n):
                if board[i][j] == 0 and is_safe(board, i, j, n):
                    board[i][j] = 1
                    queens_positions.append((i,j))
                    if solve_util(board, queens_left - 1):
                        return True
                    board[i][j] = 0
                    queens_positions.pop()
        return False
    
    if solve_util(board, 6):  # 6 remaining queens to place
        # Sort positions by row for consistent output
        queens_positions.sort()
        result = ", ".join(f"{pos[0]} {pos[1]}" for pos in queens_positions)
        print(result)
    else:
        print("No solution exists")

solve_n_queens()