def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        row = f"{i:2d} "
        for j in range(11):
            row += board[i-1][j] + " "
        print(row)

def get_liberties(x, y, board, checked):
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0
    if (x, y) in checked:
        return 0
    if board[y][x] == '.':
        return 1
    return 0

def analyze_position():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],  # 1
        ['.', '.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],  # 2
        ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.'],  # 3
        ['.', '.', '.', '.', '.', '.', 'O', 'X', 'X', '.', '.'],  # 4
        ['O', '.', '.', '.', '.', '.', '.', 'O', '.', '.', 'O'],  # 5
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 6
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 7
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
        ['.', '.', '.', '.', 'O', '.', '.', '.', 'O', '.', '.'],  # 9
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 10
        ['.', '.', '.', 'X', '.', '.', 'O', '.', 'O', '.', '.'],  # 11
    ]
    
    # Check liberties of white stones at G4, H3, I3, H2
    critical_points = [(6,3), (7,2), (8,2), (7,1)]  # G4, H3, I3, H2
    
    for x, y in critical_points:
        liberties = get_liberties(x, y-1, board, set()) + \
                   get_liberties(x+1, y, board, set()) + \
                   get_liberties(x, y+1, board, set()) + \
                   get_liberties(x-1, y, board, set())
        print(f"Liberties at {chr(65+x)}{y+1}: {liberties}")
    
    # Key point analysis
    key_points = [(6,2), (7,3), (8,3)]  # G3, H4, I4
    print("\nPotential key points:")
    for x, y in key_points:
        if board[y][x] == '.':
            print(f"{chr(65+x)}{y+1}")

analyze_position()