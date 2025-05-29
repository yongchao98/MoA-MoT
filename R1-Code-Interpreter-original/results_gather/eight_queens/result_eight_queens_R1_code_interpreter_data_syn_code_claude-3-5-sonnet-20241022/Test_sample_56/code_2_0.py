def is_safe(board, row, col, n):
    # Check row on left side
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(n):
        for j in range(n):
            # Check if position (i,j) is in same diagonal
            if (i + j == row + col) or (i - j == row - col):
                if board[i][j] == 1:
                    return False
    
    return True

def solve_n_queens(n=8):
    board = [[0 for x in range(n)] for y in range(n)]
    # Place the pre-existing queen
    board[7][5] = 1
    queens = [(7,5)]
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
        
        # Try all rows for each column
        for col in range(n):
            for row in range(n):
                # Skip if this is the pre-placed queen's position
                if row == 7 and col == 5:
                    continue
                    
                # Check if we can place a queen here
                if board[row][col] == 0:  # Position is empty
                    # Temporarily place queen
                    old_board = [row[:] for row in board]
                    board[row][col] = 1
                    
                    if is_safe(board, row, col, n):
                        queens.append((row, col))
                        if solve_util(board, queens_placed + 1):
                            return True
                        queens.pop()
                    
                    # Backtrack
                    board = old_board
        
        return False
    
    if solve_util(board, 1):  # Start with 1 as we already have one queen
        queens.sort()  # Sort by row number
        result = ", ".join(f"{row} {col}" for row, col in queens)
        print(f"<<{result}>>")
    else:
        print("No solution exists")

solve_n_queens()