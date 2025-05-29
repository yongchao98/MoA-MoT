def analyze_position(board):
    def is_valid(x, y):
        return 0 <= x < 13 and 0 <= y < 13
    
    def get_group_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        if not is_valid(x, y) or board[y][x] != color or (x,y) in visited:
            return set()
        
        visited.add((x,y))
        liberties = set()
        
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            nx, ny = x + dx, y + dy
            if is_valid(nx, ny):
                if board[ny][nx] == '.':
                    liberties.add((nx,ny))
                elif board[ny][nx] == color:
                    liberties.update(get_group_liberties(nx, ny, color, visited))
        
        return liberties

    # Create a copy of the board and analyze the white group's weaknesses
    white_group_pos = [(1,1), (1,2), (1,3), (2,3), (0,3)]  # B1, B2, B3, C3, A3
    
    # Get current liberties of the white group
    all_liberties = set()
    for x, y in white_group_pos:
        all_liberties.update(get_group_liberties(x, y, 'O'))
    
    print("Current white group liberties:", len(all_liberties))
    print("Liberty points:", [(chr(65+x), y+1) for x,y in all_liberties])
    
    # Test critical points
    critical_points = [(0,2), (2,1), (2,0)]  # A3, C2, C1
    for x, y in critical_points:
        board_copy = [row[:] for row in board]
        if board_copy[y][x] == '.':
            board_copy[y][x] = 'X'
            new_liberties = set()
            for wx, wy in white_group_pos:
                new_liberties.update(get_group_liberties(wx, wy, 'O'))
            print(f"After Black plays {chr(65+x)}{y+1}:")
            print(f"Remaining liberties: {len(new_liberties)}")
            print(f"Liberty points: {[(chr(65+x), y+1) for x,y in new_liberties]}")

# Create the board
board = [['.'] * 13 for _ in range(13)]

# Place the stones
black_stones = [(1,12), (0,4), (2,4), (3,3), (5,3), (2,2)]
white_stones = [(1,1), (1,2), (1,3), (2,3), (0,3)]

for x, y in black_stones:
    board[y][x] = 'X'
for x, y in white_stones:
    board[y][x] = 'O'

analyze_position(board)