def check_move_impact(board, move_x, move_y):
    # Check if the move at (move_x, move_y) creates a capturing threat
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    white_neighbors = 0
    black_support = 0
    
    for dx, dy in directions:
        new_x, new_y = move_x + dx, move_y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10:
            if board[new_x][new_y] == 'O':
                white_neighbors += 1
            elif board[new_x][new_y] == 'X':
                black_support += 1
    
    return white_neighbors, black_support

# Initialize the board
board = [
    ['.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.'],  # 10
    ['X', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],  # 9
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 7
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],  # 6
    ['.', 'O', '.', '.', '.', '.', 'O', 'X', 'X', '.'],  # 5
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 4
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 2
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.']   # 1
]

# Test potential moves
critical_points = [
    (4, 2, "C6"),  # C6
    (4, 4, "E6"),  # E6
    (5, 3, "D5")   # D5
]

for x, y, move_name in critical_points:
    white_neighbors, black_support = check_move_impact(board, x, y)
    print(f"\nMove at {move_name}:")
    print(f"Adjacent white stones: {white_neighbors}")
    print(f"Adjacent black support: {black_support}")