def create_board():
    board = [['.'] * 12 for _ in range(12)]
    # Place black stones (X)
    black_stones = [(11,1), (10,6), (10,9), (8,9), (7,7), (7,9), (6,6), (6,10), 
                    (5,7), (5,9), (4,8), (3,11), (2,7)]
    # Place white stones (O)
    white_stones = [(9,7), (7,6), (6,7), (6,8), (6,9), (5,8), (4,1), (2,1), (2,3)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    
    return board

def find_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return set()
    
    if board[x][y] != color:
        if board[x][y] == '.':
            return {(x,y)}
        return set()
    
    if (x,y) in visited:
        return set()
    
    visited.add((x,y))
    liberties = set()
    
    # Check all adjacent points
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_liberties = find_liberties(board, x+dx, y+dy, color, visited)
        liberties.update(new_liberties)
    
    return liberties

def analyze_position():
    board = create_board()
    # Find liberties of white group around H6
    white_liberties = find_liberties(board, 6, 7, 'O')
    print("White group liberties:", white_liberties)
    
    # Check key points around the white group
    key_points = [(5,6), (7,5), (7,8)]
    print("\nAnalyzing key points:")
    for x, y in key_points:
        if board[x][y] == '.':
            print(f"Point ({x},{y}) is empty")

analyze_position()