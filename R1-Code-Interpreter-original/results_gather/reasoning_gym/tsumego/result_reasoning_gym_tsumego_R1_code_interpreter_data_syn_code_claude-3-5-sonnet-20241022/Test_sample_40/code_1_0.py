def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = f" {9-i}"
        for j in range(9):
            row += f" {board[i][j]}"
        print(row)

def check_capture(move):
    # Initial board position (corrected)
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 9
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
        ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],  # 7
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 6
        ['.', '.', '.', '.', '.', 'X', 'O', '.', '.'],  # 5
        ['.', '.', '.', 'X', '.', '.', '.', 'O', '.'],  # 4
        ['O', '.', 'X', 'O', '.', 'O', '.', '.', '.'],  # 3
        ['.', '.', 'X', 'O', 'O', 'X', '.', '.', 'X'],  # 2
        ['.', '.', 'X', 'X', 'X', '.', 'O', '.', '.'],  # 1
    ]
    
    # Make the move
    row, col = move
    board[row][col] = 'X'
    
    print("\nAfter Black's move at E7:")
    print_board(board)
    
    # Analyze liberties of white stones at D8, E8, F7
    white_group = [(2,3), (2,4), (3,5)]  # D8, E8, F7
    liberties = set()
    for w_row, w_col in white_group:
        for dr, dc in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_row, new_col = w_row + dr, w_col + dc
            if 0 <= new_row < 9 and 0 <= new_col < 9:
                if board[new_row][new_col] == '.':
                    liberties.add((new_row, new_col))
    
    print(f"\nLiberties remaining: {len(liberties)}")
    return len(liberties) == 0

# Test the key move at E7
move = (2, 4)  # E7 in zero-based coordinates
result = check_capture(move)