def find_group(board, x, y, color, group=None, visited=None):
    if group is None:
        group = set()
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return
    
    pos = (x, y)
    if pos in visited or board[x][y] != color:
        return
        
    visited.add(pos)
    group.add(pos)
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        find_group(board, x+dx, y+dy, color, group, visited)
    return group

def find_liberties(board, group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
                liberties.add((nx, ny))
    return liberties

# Create board representation
board = [['.'] * 12 for _ in range(12)]
position = """
   A B C D E F G H I J K L
12 . X O X . . O . . . . .
11 X O O . X . . . . . . X
10 . X O X X . . . . . . .
 9 . . X . . . . . . . . .
 8 . . . . . . . . . . . .
 7 . . . . . . . . . . . .
 6 . . . . . . . . . . . .
 5 . . . . . . . . . . . .
 4 . . . . . . . . . . . .
 3 . . . . . . X . . . X .
 2 . . . X . . . O . . . O
 1 . . . . . . . . O . X .
"""

# Parse the position into the board
lines = position.strip().split('\n')
for i, line in enumerate(lines[1:], 0):
    if len(line) < 4:
        continue
    row = line[3:].split()
    for j, stone in enumerate(row):
        if stone in ['X', 'O']:
            board[i][j] = stone

# Find all white groups and their liberties
for i in range(12):
    for j in range(12):
        if board[i][j] == 'O':
            group = find_group(board, i, j, 'O')
            if group:
                liberties = find_liberties(board, group)
                if len(liberties) == 1:  # If group has only one liberty
                    liberty = list(liberties)[0]
                    print(f"Capturing move at column {chr(65+liberty[1])}, row {12-liberty[0]}")