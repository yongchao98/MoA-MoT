def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = f" {9-i}"
        for j in range(9):
            row += f" {board[i][j]}"
        print(row)

def check_capture(move):
    # Initial board position
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', 'X', '.', 'O', '.', '.', 'O', '.'],
        ['.', '.', '.', '.', '.', 'X', 'O', '.', '.'],
        ['.', '.', '.', 'X', '.', '.', '.', 'O', '.'],
        ['O', '.', 'X', 'O', '.', 'O', '.', '.', '.'],
        ['.', '.', 'X', 'O', 'O', 'X', '.', '.', 'X'],
        ['.', '.', 'X', 'X', 'X', '.', 'O', '.', '.']
    ]
    
    # Make the move
    row, col = move
    board[row][col] = 'X'
    
    print("\nAfter Black's move at E7:")
    print_board(board)
    
    # Check if white stones are captured (simplified check)
    white_stones = [(7,3), (7,4), (6,5)]  # D8, E8, F7
    surrounded = True
    for w_row, w_col in white_stones:
        for dr, dc in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_row, new_col = w_row + dr, w_col + dc
            if 0 <= new_row < 9 and 0 <= new_col < 9:
                if board[new_row][new_col] == '.':
                    surrounded = False
                    break
    
    return surrounded

# Test the key move at E7
move = (6, 4)  # E7 in zero-based coordinates
result = check_capture(move)
print(f"\nCan capture white stones: {result}")