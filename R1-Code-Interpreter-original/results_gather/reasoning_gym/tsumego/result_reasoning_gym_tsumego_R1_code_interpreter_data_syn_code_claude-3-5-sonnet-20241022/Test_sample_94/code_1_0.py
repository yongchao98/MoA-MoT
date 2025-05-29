def create_board():
    board = [['.'] * 13 for _ in range(13)]
    # White stones (O)
    white_positions = [(1,6), (1,5), (6,6), (1,7), (8,6), (9,12)]
    # Black stones (X)
    black_positions = [(0,7), (1,8), (2,7), (2,5), (0,5), (2,4), (3,6), 
                      (6,6), (10,1), (10,2), (10,10), (10,11), (4,11), (12,7)]
    
    for x, y in white_positions:
        board[y][x] = 'O'
    for x, y in black_positions:
        board[y][x] = 'X'
    return board

def get_group_liberties(board, start_x, start_y, visited=None):
    if visited is None:
        visited = set()
    
    color = board[start_y][start_x]
    if color == '.':
        return set()
    
    liberties = set()
    stack = [(start_x, start_y)]
    stone_group = set()
    
    while stack:
        x, y = stack.pop()
        if (x, y) in visited:
            continue
        
        visited.add((x, y))
        stone_group.add((x, y))
        
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 13 and 0 <= new_y < 13:
                if board[new_y][new_x] == '.':
                    liberties.add((new_x, new_y))
                elif board[new_y][new_x] == color and (new_x, new_y) not in visited:
                    stack.append((new_x, new_y))
    
    return liberties, stone_group

def analyze_groups():
    board = create_board()
    visited = set()
    
    for y in range(13):
        for x in range(13):
            if board[y][x] == 'O' and (x,y) not in visited:
                liberties, group = get_group_liberties(board, x, y, visited)
                if len(liberties) <= 2:  # Check groups with 2 or fewer liberties
                    print(f"White group at {[(x+1,13-y) for x,y in group]}")
                    print(f"Has liberties at: {[(chr(65+x),13-y) for x,y in liberties]}")

analyze_groups()