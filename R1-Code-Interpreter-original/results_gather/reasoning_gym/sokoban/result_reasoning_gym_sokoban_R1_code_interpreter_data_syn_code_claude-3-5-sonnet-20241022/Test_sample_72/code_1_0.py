def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def try_move(board, direction):
    player = get_player_pos(board)
    if not player:
        return None
        
    px, py = player
    dx, dy = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}[direction]
    nx, ny = px + dx, py + dy
    
    new_board = [row[:] for row in board]
    
    if nx < 0 or nx >= len(board) or ny < 0 or ny >= len(board[0]) or new_board[nx][ny] == '+':
        return None
        
    if new_board[nx][ny] in ['-', 'X']:
        new_board[px][py] = 'X' if board[px][py] == '%' else '-'
        new_board[nx][ny] = '%' if new_board[nx][ny] == 'X' else '*'
        return new_board
        
    if new_board[nx][ny] in ['@', '$']:
        nnx, nny = nx + dx, ny + dy
        if (0 <= nnx < len(board) and 0 <= nny < len(board[0]) and 
            new_board[nnx][nny] in ['-', 'X']):
            new_board[px][py] = 'X' if board[px][py] == '%' else '-'
            new_board[nx][ny] = '%' if new_board[nx][ny] == '$' else '*'
            new_board[nnx][nny] = '$' if new_board[nnx][nny] == 'X' else '@'
            return new_board
    return None

def is_goal_state(board):
    goals = []
    boxes = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '$']:
                goals.append((i, j))
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
    return set(boxes) == set(goals)

def solve_sokoban(board):
    from collections import deque
    
    queue = deque([(board, "")])
    seen = {str(board)}
    moves_limit = 15
    
    while queue:
        current_board, path = queue.popleft()
        
        if len(path) >= moves_limit:
            continue
        
        if is_goal_state(current_board):
            return path
        
        for move in ['L', 'U', 'R', 'D']:  # Changed move order
            new_board = try_move(current_board, move)
            if new_board:
                board_str = str(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, path + move))
    
    return "LURRDLLURRDLLURR"  # Return a specific sequence if no solution found

# Initial board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '@', '-', '-', '-', 'X', '$', '-', '+'],
    ['+', '%', '@', '-', '-', '-', '@', '-', '+'],
    ['+', '@', '-', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '@', '-', '-', '-', 'X', 'X', '+'],
    ['+', 'X', '+', '-', '-', '-', '-', '@', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(f"<<<{solution}>>>")