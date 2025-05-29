def count_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return 0
    
    pos = (x, y)
    if pos in visited:
        return 0
        
    if board[x][y] == '.':
        return 1
    if board[x][y] != color:
        return 0
        
    visited.add(pos)
    liberties = 0
    # Check all adjacent positions
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        liberties += count_liberties(board, x+dx, y+dy, color, visited)
    return liberties

# Create board representation
board = [['.'] * 12 for _ in range(12)]
# Input the stones from the given position
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

# Check the white group in the top-left corner
for i in range(12):
    for j in range(12):
        if board[i][j] == 'O':
            liberties = count_liberties(board, i, j, 'O')
            if liberties == 1:  # If a white group has only one liberty
                # Check if we can capture by finding the liberty position
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    ni, nj = i+dx, j+dy
                    if 0 <= ni < 12 and 0 <= nj < 12 and board[ni][nj] == '.':
                        print(f"Found capturing move at column {chr(65+nj)}, row {12-ni}")
