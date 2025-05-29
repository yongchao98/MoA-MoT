def analyze_move(board, move_x, move_y):
    def is_valid(x, y):
        return 0 <= x < 9 and 0 <= y < 9

    def get_group_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        
        if not is_valid(x, y) or (x, y) in visited:
            return set()
        
        if board[y][x] != color:
            return {(x, y)} if board[y][x] == '.' else set()
            
        visited.add((x, y))
        liberties = set()
        
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if is_valid(new_x, new_y):
                liberties.update(get_group_liberties(new_x, new_y, color, visited))
        
        return liberties

    # Make move
    test_board = [row[:] for row in board]
    test_board[move_y][move_x] = 'X'
    
    # Check upper white group (G7)
    upper_liberties = get_group_liberties(6, 2, 'O')  # G7
    
    # Check lower white group (F4)
    lower_liberties = get_group_liberties(5, 5, 'O')  # F4
    
    return {
        'move': f"{chr(65+move_x)}{9-move_y}",
        'upper_liberties': len(upper_liberties),
        'upper_liberty_points': sorted(list((chr(65+x), 9-y) for x, y in upper_liberties)),
        'lower_liberties': len(lower_liberties),
        'lower_liberty_points': sorted(list((chr(65+x), 9-y) for x, y in lower_liberties))
    }

# Initialize board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'O', 'X'],
    ['X', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', 'X', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.']
]

# Analyze key moves
moves_to_check = [(6,4), (6,3), (5,4)]  # G5, G6, F5
for x, y in moves_to_check:
    result = analyze_move(board, x, y)
    print(f"\nIf Black plays {result['move']}:")
    print(f"Upper white group liberties: {result['upper_liberties']}")
    print(f"Upper liberty points: {result['upper_liberty_points']}")
    print(f"Lower white group liberties: {result['lower_liberties']}")
    print(f"Lower liberty points: {result['lower_liberty_points']}")