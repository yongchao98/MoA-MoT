def check_capture_sequence():
    # Initialize board
    board = [
        ['O', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],  # 10
        ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],  # 9
        ['.', 'O', '.', '.', 'O', '.', '.', '.', '.', '.'],  # 8
        ['X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 7
        ['O', 'O', '.', 'O', '.', 'X', '.', '.', '.', '.'],  # 6
        ['X', 'O', 'X', 'O', '.', '.', '.', '.', '.', '.'],  # 5
        ['.', 'X', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 4
        ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],  # 3
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', '.'],  # 2
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1
    ]
    
    def get_group_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
            
        if (x, y) in visited or x < 0 or x >= 10 or y < 0 or y >= 10:
            return set()
            
        if board[y][x] != color:
            if board[y][x] == '.':
                return {(x, y)}
            return set()
            
        visited.add((x, y))
        liberties = set()
        for dx, dy in [(1,0), (-1,0), (0,1), (0,-1)]:
            liberties.update(get_group_liberties(x+dx, y+dy, color, visited))
        return liberties

    # Test C6 move
    x, y = 2, 6  # C6
    board[y][x] = 'X'
    
    # Check liberties of adjacent white groups
    white_groups = [(1,4), (1,6)]  # B5, B6
    for wx, wy in white_groups:
        liberties = get_group_liberties(wx, wy, 'O')
        print(f"White group at {chr(65+wx)}{10-wy} has {len(liberties)} liberties after Black C6")
        print(f"Liberties at: {[(chr(65+lx), 10-ly) for lx, ly in liberties]}")

check_capture_sequence()