def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Place stones exactly as shown
    black = "A8 B9 B7 C10 C6 D9 D7 H9 J7 J6 J5 J4 B1 I1"
    white = "B8 C8 C7 D3 F5 I5 I2"
    
    for pos in black.split():
        col = ord(pos[0]) - ord('A')
        row = 10 - int(pos[1:])
        board[row][col] = 'X'
    for pos in white.split():
        col = ord(pos[0]) - ord('A')
        row = 10 - int(pos[1:])
        board[row][col] = 'O'
    return board

def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if (x, y) in visited or x < 0 or x >= 10 or y < 0 or y >= 10:
        return set()
    
    if board[x][y] == '.':
        return {(x, y)}
    
    if board[x][y] != 'O':  # We only care about white stones
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10:
            if board[nx][ny] == '.':
                liberties.add((nx, ny))
            elif board[nx][ny] == 'O':
                liberties.update(count_liberties(board, nx, ny, visited))
    
    return liberties

def analyze_white_groups():
    board = create_board()
    white_groups = [
        # Upper group
        [(1,7), (2,7), (2,6)],  # B8,C8,C7 in board coordinates
        # Other stones
        [(3,3)],  # D3
        [(5,5)],  # F5
        [(5,8)],  # I5
        [(2,8)]   # I2
    ]
    
    for group in white_groups:
        liberties = set()
        stones = set()
        for x, y in group:
            stones.update({(x,y)})
            liberties.update(count_liberties(board, x, y, set()))
        
        # Convert to board notation
        stones_notation = [f"{chr(y + ord('A'))}{10-x}" for x, y in stones]
        liberties_notation = [f"{chr(y + ord('A'))}{10-x}" for x, y in liberties]
        print(f"Group {stones_notation} has liberties at {liberties_notation}")

analyze_white_groups()