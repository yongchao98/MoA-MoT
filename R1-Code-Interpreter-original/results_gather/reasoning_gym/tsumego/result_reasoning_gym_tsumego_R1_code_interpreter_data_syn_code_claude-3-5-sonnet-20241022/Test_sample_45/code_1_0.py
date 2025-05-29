def print_board_state():
    # Create board (11x11)
    board = [
        ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 11 (0)
        ['.', 'O', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 10 (1)
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 9  (2)
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', 'X', 'X'],  # 8  (3)
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],  # 7  (4)
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', 'O', '.'],  # 6  (5)
        ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 5  (6)
        ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],  # 4  (7)
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O'],  # 3  (8)
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 2  (9)
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1  (10)
    ]
    
    # Check I7 position and surroundings
    i7_row = 10 - 7  # Convert from Go coordinates to array index
    i7_col = ord('I') - ord('A')
    
    print(f"I7 position [row={i7_row}, col={i7_col}]: {board[i7_row][i7_col]}")
    print("\nSurrounding stones:")
    print(f"J7 [row={i7_row}, col={i7_col+1}]: {board[i7_row][i7_col+1]}")
    print(f"K7 [row={i7_row}, col={i7_col+2}]: {board[i7_row][i7_col+2]}")
    print(f"J8 [row={i7_row-1}, col={i7_col+1}]: {board[i7_row-1][i7_col+1]}")
    print(f"J6 [row={i7_row+1}, col={i7_col+1}]: {board[i7_row+1][i7_col+1]}")

print_board_state()