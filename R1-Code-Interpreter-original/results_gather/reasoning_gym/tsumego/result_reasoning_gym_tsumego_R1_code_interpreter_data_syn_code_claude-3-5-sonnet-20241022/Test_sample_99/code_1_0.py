def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0
    
    if board[x][y] == '.':
        return 1
    
    if (x, y) in visited:
        return 0
        
    if board[x][y] == 'O':
        visited.add((x, y))
        liberties = 0
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            liberties += count_liberties(board, x+dx, y+dy, visited)
        return liberties
    
    return 0

def would_capture(board, move_x, move_y):
    if board[move_x][move_y] != '.':
        return 0
        
    # Make the move
    board[move_x][move_y] = 'X'
    captured = 0
    
    # Check all adjacent positions
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        x, y = move_x + dx, move_y + dy
        if 0 <= x < 11 and 0 <= y < 11 and board[x][y] == 'O':
            if count_liberties(board, x, y) == 0:
                # Count connected white stones that would be captured
                temp_captured = 1
                stack = [(x, y)]
                visited = {(x, y)}
                while stack:
                    cx, cy = stack.pop()
                    for ndx, ndy in [(0,1), (0,-1), (1,0), (-1,0)]:
                        nx, ny = cx + ndx, cy + ndy
                        if (0 <= nx < 11 and 0 <= ny < 11 and 
                            board[nx][ny] == 'O' and (nx, ny) not in visited):
                            temp_captured += 1
                            visited.add((nx, ny))
                            stack.append((nx, ny))
                captured += temp_captured
    
    # Undo the move
    board[move_x][move_y] = '.'
    return captured

# Initialize the board
board = [['.'] * 11 for _ in range(11)]
# Set up the board according to the given position
positions = {
    'X': [(10,3), (9,2), (9,4), (8,2), (8,4), (8,8), (8,10), (7,0), (7,3), (6,8), (4,4)],
    'O': [(9,0), (9,2), (9,3), (8,3), (9,5), (9,7), (4,7), (3,10)]
}

for piece, coords in positions.items():
    for x, y in coords:
        board[x][y] = piece

# Try all possible moves and find the one that captures the most stones
max_capture = 0
best_move = None

for i in range(11):
    for j in range(11):
        if board[i][j] == '.':
            captures = would_capture(board, i, j)
            if captures > max_capture:
                max_capture = captures
                best_move = (j, 10-i)  # Convert to board coordinates

if best_move:
    print(f"Best move is {chr(65+best_move[0])}{best_move[1]+1} capturing {max_capture} stones")