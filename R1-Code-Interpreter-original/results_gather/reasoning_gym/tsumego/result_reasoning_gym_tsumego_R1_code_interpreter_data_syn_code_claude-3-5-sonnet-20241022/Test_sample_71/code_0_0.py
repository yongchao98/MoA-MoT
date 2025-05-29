def print_board(board):
    print("   A B C D E F G H I J K L")
    for i in range(12):
        row = 12 - i
        print(f"{row:2d}", end=" ")
        for j in range(12):
            print(board[i][j], end=" ")
        print()

def analyze_position():
    board = [['.'] * 12 for _ in range(12)]
    
    # Place black stones (X)
    black_stones = [(0,0), (0,11), (2,11), (7,4), (8,10), (9,2), (9,3), (10,1), (10,6), (10,8), (11,9)]
    
    # Place white stones (O)
    white_stones = [(2,3), (5,7), (5,11), (9,2), (9,3), (10,2), (10,3), (11,2), (11,3)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    
    # Analyze vital points
    vital_points = [(9,1), (9,2), (8,2), (8,3)]  # Potential key moves
    
    print("Original position:")
    print_board(board)
    print("\nPotential vital points (coordinates in board matrix):")
    for point in vital_points:
        print(f"({point[0]}, {point[1]})")

analyze_position()