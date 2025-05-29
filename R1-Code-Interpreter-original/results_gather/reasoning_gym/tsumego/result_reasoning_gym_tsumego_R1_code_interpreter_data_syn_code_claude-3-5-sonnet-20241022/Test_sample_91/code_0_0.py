def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return 0, set()
    
    if board[y][x] == '.':
        return 1, {(x,y)}
    
    if (x,y) in visited:
        return 0, set()
        
    if board[y][x] != 'O':  # if not white
        return 0, set()
        
    visited.add((x,y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        count, libs = count_liberties(board, new_x, new_y, visited)
        liberties.update(libs)
    
    return len(liberties), liberties

# Create the board
board = [
    ['.', '.', '.', '.', 'O', '.', 'O', '.', '.'],  # 9
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', 'X', '.', '.', '.', '.', 'O', '.'],  # 7
    ['.', '.', '.', '.', '.', '.', 'X', '.', 'X'],  # 6
    ['.', 'O', '.', '.', '.', 'X', 'O', 'O', 'O'],  # 5
    ['.', '.', '.', 'X', '.', 'O', 'X', 'O', 'X'],  # 4
    ['.', 'O', '.', '.', '.', '.', '.', 'X', '.'],  # 3
    ['.', 'O', '.', 'O', 'X', '.', '.', '.', '.'],  # 2
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1
]

# Check each white stone group and its liberties
best_move = None
max_group_size = 0

for y in range(9):
    for x in range(9):
        if board[y][x] == 'O':
            visited = set()
            liberties_count, liberty_positions = count_liberties(board, x, y)
            if liberties_count == 1 and len(visited) > max_group_size:
                max_group_size = len(visited)
                best_move = list(liberty_positions)[0]

if best_move:
    col = chr(ord('A') + best_move[0])
    row = 9 - best_move[1]
    print(f"{col}{row} can capture {max_group_size} stones")