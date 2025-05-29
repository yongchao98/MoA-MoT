def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0
    
    if board[x][y] == '.':
        return 1
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    liberties = 0
    
    # Check all adjacent positions
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 11 and 0 <= new_y < 11:
            if board[new_x][new_y] == '.':
                liberties += 1
            elif board[new_x][new_y] == board[x][y]:
                liberties += count_liberties(board, new_x, new_y, visited)
    
    return liberties

# Create the board
board = [['.'] * 11 for _ in range(11)]

# Place the stones
stones = {
    'X': [(10,4), (9,4), (8,9), (7,3), (6,3), (5,8), (4,2), (2,9), (2,8), (2,7)],
    'O': [(9,1), (7,2), (6,9), (6,10), (2,10)]
}

for color in stones:
    for x, y in stones[color]:
        board[x][y] = color

# Find the most vulnerable white group
min_liberties = float('inf')
vulnerable_group = None

for x in range(11):
    for y in range(11):
        if board[x][y] == 'O':
            liberties = count_liberties(board, x, y)
            if liberties < min_liberties:
                min_liberties = liberties
                vulnerable_group = (x, y)

# The move should be adjacent to the most vulnerable group
if vulnerable_group:
    x, y = vulnerable_group
    print(f"Most vulnerable white group at: {chr(65+y)}{11-x}")
    print(f"Best capturing move: {chr(65+y-1)}{11-x}")