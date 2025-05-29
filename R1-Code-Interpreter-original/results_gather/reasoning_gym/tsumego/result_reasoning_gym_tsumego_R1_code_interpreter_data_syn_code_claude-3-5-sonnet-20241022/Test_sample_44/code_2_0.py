def create_board():
    board = [['.'] * 11 for _ in range(11)]
    
    # Place Black stones (X)
    black_positions = [
        (4,10), (6,10),  # row 10
        (3,9), (7,9),    # row 9
        (4,8), (6,8),    # row 8
        (0,7),           # row 7
        (4,6),           # row 6
        (3,5),           # row 5
        (0,4), (6,4), (7,4),  # row 4
        (0,3),           # row 3
        (6,1)            # row 1
    ]
    
    # Place White stones (O)
    white_positions = [
        (5,11),          # row 11
        (4,9), (5,9), (6,9),  # row 9
        (5,8), (10,8),   # row 8
        (9,5)            # row 5
    ]
    
    for x, y in black_positions:
        board[11-y][x] = 'X'
    for x, y in white_positions:
        board[11-y][x] = 'O'
    
    return board

def find_connected_group(board, start_x, start_y):
    if board[start_y][start_x] != 'O':
        return set()
    
    group = set()
    liberties = set()
    stack = [(start_x, start_y)]
    
    while stack:
        x, y = stack.pop()
        if (x, y) in group:
            continue
            
        if board[y][x] == 'O':
            group.add((x, y))
            # Check neighbors
            for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                new_x, new_y = x + dx, y + dy
                if 0 <= new_x < 11 and 0 <= new_y < 11:
                    if board[new_y][new_x] == 'O' and (new_x, new_y) not in group:
                        stack.append((new_x, new_y))
                    elif board[new_y][new_x] == '.':
                        liberties.add((new_x, new_y))
    
    return group, liberties

def analyze_white_groups(board):
    analyzed = set()
    groups = []
    
    for y in range(11):
        for x in range(11):
            if board[y][x] == 'O' and (x, y) not in analyzed:
                group, liberties = find_connected_group(board, x, y)
                analyzed.update(group)
                groups.append((group, liberties))
    
    return groups

board = create_board()
groups = analyze_white_groups(board)

print("White groups analysis:")
for group, liberties in groups:
    print(f"\nGroup size: {len(group)}")
    print(f"Stones: {[(x+1, 11-y) for x,y in group]}")
    print(f"Liberties: {[(x+1, 11-y) for x,y in liberties]}")
    print(f"Number of liberties: {len(liberties)}")