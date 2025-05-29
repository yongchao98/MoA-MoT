def create_board():
    board = [['.'] * 11 for _ in range(11)]
    # Add black stones (X)
    black_stones = [(6,4), (8,4), (9,4), (7,5), (5,6), (7,7), (8,8), (9,8), (10,8), (1,7), (4,7), (3,8), (1,9)]
    # Add white stones (O)
    white_stones = [(11,6), (9,6), (9,7), (8,6), (8,8), (5,11), (5,9)]
    
    for x, y in black_stones:
        board[y-1][x-1] = 'X'
    for x, y in white_stones:
        board[y-1][x-1] = 'O'
    return board

def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0, set()
    
    if board[y][x] == '.':
        return 1, {(x,y)}
    
    if (x,y) in visited:
        return 0, set()
        
    if board[y][x] != 'O':
        return 0, set()
    
    visited.add((x,y))
    liberties = set()
    
    # Check all four directions
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        count, libs = count_liberties(board, x+dx, y+dy, visited)
        liberties.update(libs)
    
    return len(liberties), liberties

def find_vulnerable_groups(board):
    vulnerable = []
    for y in range(11):
        for x in range(11):
            if board[y][x] == 'O':
                liberties, liberty_points = count_liberties(board, x, y)
                if liberties <= 2:  # Groups with 2 or fewer liberties are vulnerable
                    vulnerable.append((x, y, liberties, liberty_points))
    return vulnerable

board = create_board()
vulnerable = find_vulnerable_groups(board)
print("Vulnerable white groups:")
for x, y, libs, points in vulnerable:
    print(f"Stone at ({x+1},{y+1}) has {libs} liberties at positions: {[(p[0]+1,p[1]+1) for p in points]}")