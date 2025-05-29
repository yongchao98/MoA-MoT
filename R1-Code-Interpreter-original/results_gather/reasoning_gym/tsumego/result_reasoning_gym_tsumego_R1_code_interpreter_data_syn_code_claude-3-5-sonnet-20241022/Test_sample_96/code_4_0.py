def analyze_capturing_potential(board, move_x, move_y, move_name):
    # Initialize the board with the test move
    test_board = [row[:] for row in board]
    test_board[move_x][move_y] = 'X'
    
    # Check if move creates immediate pressure on white group
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    pressure_points = 0
    connection_to_black = 0
    cutting_potential = 0
    
    # Check direct adjacencies
    for dx, dy in directions:
        new_x, new_y = move_x + dx, move_y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10:
            if test_board[new_x][new_y] == 'O':
                # Check if this white stone is part of the main group
                if new_x <= 4 and new_y >= 2:  # Upper right quadrant where main white group is
                    pressure_points += 2
                else:
                    pressure_points += 1
            elif test_board[new_x][new_y] == 'X':
                connection_to_black += 1
    
    # Check if move cuts white's connection paths
    for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        new_x, new_y = move_x + dx, move_y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10:
            if test_board[new_x][new_y] == 'O':
                cutting_potential += 1

    print(f"\nDetailed analysis for {move_name}:")
    print(f"Pressure on white group: {pressure_points}")
    print(f"Connections to black: {connection_to_black}")
    print(f"Cutting potential: {cutting_potential}")
    print(f"Total effectiveness score: {pressure_points + connection_to_black + cutting_potential}")

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

# Analyze critical moves
moves_to_analyze = [
    (4, 2, "C6"),  # C6
    (5, 3, "D5"),  # D5
    (4, 4, "E6")   # E6
]

for x, y, move in moves_to_analyze:
    analyze_capturing_potential(board, x, y, move)