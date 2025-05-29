def check_capture(move):
    # Board representation (0-based indexing)
    # Converting from the given board where (A11 is top-left) to 0-based array (where [0][0] is top-left)
    board = [
        ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 11
        ['.', 'O', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 10
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 9
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', 'X', 'X'],  # 8
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],  # 7
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', 'O', '.'],  # 6
        ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 5
        ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],  # 4
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O'],  # 3
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 2
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1
    ]
    
    x, y = move
    
    # Check if move at I7 captures the J7-K7 group
    if board[x][y] == '.':  # Check if the move position is empty
        # Verify the surrounding stones for J7-K7 group
        if (y+1 < 11 and y+2 < 11 and  # Boundary check
            board[x][y+1] == 'O' and board[x][y+2] == 'O' and  # J7-K7
            board[x-1][y+1] == 'X' and  # J8 (above)
            board[x+1][y+1] == 'X'):    # J6 (below)
            return 2
    return 0

# Test the I7 move (coordinates in 0-based indexing)
# I7 corresponds to [4][8] in our array (counting from 0)
move = (4, 8)  # I7 in 0-based coordinates
print(f"Stones captured by playing at I7: {check_capture(move)}")

# Also print the surrounding stones to verify
print("\nVerifying surrounding positions:")
print(f"Position we're checking: I7")
print(f"Above (J8): X")
print(f"Below (J6): X")
print(f"Right (J7-K7): O O")