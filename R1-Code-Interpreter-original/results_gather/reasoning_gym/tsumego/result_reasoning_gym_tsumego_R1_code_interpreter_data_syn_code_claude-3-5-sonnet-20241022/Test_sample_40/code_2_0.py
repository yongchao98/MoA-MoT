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
        ['.', '.', 'X', 'X', 'X', '.', 'O', '.', '.'],  # 9
        ['.', '.', 'X', 'O', 'O', 'X', '.', '.', 'X'],  # 8
        ['O', '.', 'X', 'O', '.', 'O', '.', '.', '.'],  # 7
        ['.', '.', '.', 'X', '.', '.', '.', 'O', '.'],  # 6
        ['.', '.', '.', '.', '.', 'X', 'O', '.', '.'],  # 5
        ['.', '.', 'X', '.', 'O', '.', '.', 'O', '.'],  # 4
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 3
        ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],  # 2
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1
    ]
    
    # Make the move
    row, col = move
    board[row][col] = 'X'
    
    print("\nAfter Black's move at E7:")
    print_board(board)
    
    # Check liberties of white stones at D8, E8, F7
    white_group = [(1,3), (1,4), (2,5)]  # D8, E8, F7
    checked = set()
    liberties = set()
    
    def count_liberties(r, c):
        if (r,c) in checked or not (0 <= r < 9 and 0 <= c < 9):
            return
        if board[r][c] == '.':
            liberties.add((r,c))
            return
        if board[r][c] != 'O':
            return
        checked.add((r,c))
        for dr, dc in [(0,1), (0,-1), (1,0), (-1,0)]:
            count_liberties(r+dr, c+dc)
    
    # Start checking from each stone in the white group
    for r, c in white_group:
        count_liberties(r, c)
    
    print(f"\nLiberties remaining: {len(liberties)}")
    print(f"Liberty points: {liberties}")
    return len(liberties) == 1  # If only one liberty remains at E7

# Test the key move at E7
move = (2, 4)  # E7 in zero-based coordinates
result = check_capture(move)