def analyze_move(board, move_x, move_y):
    def is_valid(x, y):
        return 0 <= x < 13 and 0 <= y < 13
    
    def get_group_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        if not is_valid(x, y) or (x,y) in visited or board[y][x] != color:
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

    # Make a copy and play the move
    board_copy = [row[:] for row in board]
    if board_copy[move_y][move_x] != '.':
        return "Invalid move - position occupied"
    
    board_copy[move_y][move_x] = 'X'
    
    # Check impact on white groups
    white_groups = []
    visited = set()
    
    # Find all adjacent white groups
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        nx, ny = move_x + dx, move_y + dy
        if is_valid(nx, ny) and board_copy[ny][nx] == 'O' and (nx,ny) not in visited:
            group = get_group_liberties(nx, ny, 'O')
            if group:
                white_groups.append(group)
                visited.update(group)
    
    # Analyze the impact
    result = f"Move at {chr(65+move_x)}{move_y+1}:\n"
    for i, group in enumerate(white_groups):
        result += f"Adjacent white group {i+1} liberties: {len(group)}\n"
        result += f"Liberty points: {[(chr(65+x), y+1) for x,y in group]}\n"
    
    return result

# Create board
board = [['.'] * 13 for _ in range(13)]

# Place stones
black_stones = [(1,12), (0,4), (2,4), (3,3), (5,3), (2,2)]
white_stones = [(1,1), (1,2), (1,3), (2,3), (0,3)]

for x, y in black_stones:
    board[y][x] = 'X'
for x, y in white_stones:
    board[y][x] = 'O'

# Test critical points
critical_points = [(0,0), (1,0), (2,0)]  # A1, B1, C1
for x, y in critical_points:
    print(analyze_move(board, x, y))