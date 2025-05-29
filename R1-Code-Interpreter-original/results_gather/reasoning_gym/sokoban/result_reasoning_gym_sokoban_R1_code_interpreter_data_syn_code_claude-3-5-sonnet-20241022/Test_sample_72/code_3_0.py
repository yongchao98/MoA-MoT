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
    
    # Copy board
    new_board = [row[:] for row in board]
    
    # Check boundaries and walls
    if nx < 0 or nx >= len(board) or ny < 0 or ny >= len(board[0]) or new_board[nx][ny] == '+':
        return None
        
    # Moving to empty space or goal
    if new_board[nx][ny] in ['-', 'X']:
        new_board[px][py] = 'X' if board[px][py] == '%' else '-'
        new_board[nx][ny] = '%' if new_board[nx][ny] == 'X' else '*'
        return new_board
        
    # Pushing box
    if new_board[nx][ny] in ['@', '$']:
        nnx, nny = nx + dx, ny + dy
        if (0 <= nnx < len(board) and 0 <= nny < len(board[0]) and 
            new_board[nnx][nny] in ['-', 'X']):
            new_board[px][py] = 'X' if board[px][py] == '%' else '-'
            new_board[nx][ny] = '%' if new_board[nx][ny] == '$' else '*'
            new_board[nnx][nny] = '$' if new_board[nnx][nny] == 'X' else '@'
            return new_board
    return None

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(board):
    from collections import deque
    
    initial_state = (board, "")
    queue = deque([initial_state])
    seen = {board_to_string(board)}
    moves_limit = 20
    
    while queue:
        current_board, path = queue.popleft()
        
        if len(path) >= moves_limit:
            continue
            
        # Check if all boxes are on goals
        all_boxes_on_goals = True
        for row in current_board:
            if '@' in row:  # If there's a box not on a goal
                all_boxes_on_goals = False
                break
        
        if all_boxes_on_goals:
            return path
            
        # Try all possible moves
        for move in ['U', 'D', 'L', 'R']:
            new_board = try_move(current_board, move)
            if new_board:
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, path + move))
                    
    return None

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