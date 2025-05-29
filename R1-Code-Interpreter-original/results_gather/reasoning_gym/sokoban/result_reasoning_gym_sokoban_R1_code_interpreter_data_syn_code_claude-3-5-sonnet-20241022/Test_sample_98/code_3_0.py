def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def is_valid_move(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def make_move(board, player_pos, direction):
    px, py = player_pos
    nx, ny = get_next_pos(px, py, direction)
    
    if not is_valid_move(board, nx, ny):
        return None
        
    new_board = [list(row) for row in board]
    
    # Moving to empty space
    if new_board[nx][ny] == '-':
        new_board[px][py] = '-'
        new_board[nx][ny] = '*'
        return (new_board, (nx, ny))
    
    # Moving to goal
    if new_board[nx][ny] == 'X':
        new_board[px][py] = '-'
        new_board[nx][ny] = '*'
        return (new_board, (nx, ny))
    
    # Pushing box
    if new_board[nx][ny] in ['@', '$']:
        bx, by = get_next_pos(nx, ny, direction)
        if not is_valid_move(board, bx, by) or new_board[bx][by] in ['@', '$']:
            return None
            
        new_board[px][py] = '-'
        new_board[nx][ny] = '*'
        new_board[bx][by] = '$' if new_board[bx][by] == 'X' else '@'
        return (new_board, (nx, ny))
    
    return None

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_solved(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell == 'X':
                goals += 1
            elif cell == '$':
                boxes_on_goals += 1
    return goals == 0 and boxes_on_goals > 0

def solve(board, start_pos, max_depth=20):
    from collections import deque
    
    queue = deque([(board, start_pos, "")])
    visited = set()
    
    while queue:
        current_board, pos, path = queue.popleft()
        
        if len(path) > max_depth:
            continue
            
        if is_solved(current_board):
            return path
            
        state = (board_to_string(current_board), pos)
        if state in visited:
            continue
            
        visited.add(state)
        
        for direction in ['U', 'R', 'D', 'L']:
            result = make_move(current_board, pos, direction)
            if result:
                new_board, new_pos = result
                queue.append((new_board, new_pos, path + direction))
    
    return None

# Initial board
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '-', '@', 'X', '@', '-', '-', '+'],
    ['+', 'X', '-', '*', 'X', '-', '-', '+'],
    ['+', '-', '@', 'X', 'X', '@', '-', '+'],
    ['+', '@', '-', '-', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+']
]

# Find player position
player_pos = None
for i in range(len(board)):
    for j in range(len(board[0])):
        if board[i][j] == '*':
            player_pos = (i, j)
            break
    if player_pos:
        break

solution = solve(board, player_pos)
print(solution)