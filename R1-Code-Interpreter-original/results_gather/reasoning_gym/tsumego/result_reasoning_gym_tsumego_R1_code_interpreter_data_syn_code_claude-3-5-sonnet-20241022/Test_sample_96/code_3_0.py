def analyze_position(board, move_x, move_y, move_name):
    # Check surrounding area influence and potential ladder/net formations
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    extended_area = []
    
    # Check 2-space radius around the move
    for dx, dy in directions:
        for dist in [1, 2]:
            new_x, new_y = move_x + dx * dist, move_y + dy * dist
            if 0 <= new_x < 10 and 0 <= new_y < 10:
                extended_area.append(board[new_x][new_y])
    
    # Count stones in extended area
    black_count = extended_area.count('X')
    white_count = extended_area.count('O')
    empty_count = extended_area.count('.')
    
    print(f"\nAnalysis for move {move_name}:")
    print(f"Black stones in area: {black_count}")
    print(f"White stones in area: {white_count}")
    print(f"Empty points: {empty_count}")
    
    # Check if move creates a ladder or net
    creates_ladder = False
    if black_count >= 2 and white_count >= 1:
        for i in range(len(extended_area)-1):
            if extended_area[i] == 'X' and extended_area[i+1] == 'O':
                creates_ladder = True
                break
    
    if creates_ladder:
        print("This move potentially creates a ladder/net")

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

# Analyze key moves
critical_moves = [
    (4, 2, "C6"),
    (5, 3, "D5"),
    (4, 3, "E6")
]

for x, y, move in critical_moves:
    analyze_position(board, x, y, move)