from collections import deque

def get_player_pos(board):
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell in ['*', '%']:
                return (i, j)
    return None

def make_move(board, direction):
    board = [list(row) for row in board]
    player = get_player_pos(board)
    if not player:
        return None
        
    px, py = player
    dx, dy = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}[direction]
    nx, ny = px + dx, py + dy
    
    # Check if move is within bounds
    if nx < 0 or ny < 0 or nx >= len(board) or ny >= len(board[0]):
        return None
    
    # Check if hitting wall
    if board[nx][ny] == '+':
        return None
        
    # Moving to empty space or goal
    if board[nx][ny] in ['-', 'X']:
        # Move player
        board[nx][ny] = '*' if board[nx][ny] == '-' else '%'
        board[px][py] = '-' if board[px][py] == '*' else 'X'
        return board
        
    # Moving box
    if board[nx][ny] in ['@', '$']:
        bx, by = nx + dx, ny + dy
        
        # Check if box move is valid
        if (bx < 0 or by < 0 or bx >= len(board) or by >= len(board[0]) or
            board[bx][by] in ['+', '@', '$']):
            return None
            
        # Move box
        board[bx][by] = '@' if board[bx][by] == '-' else '$'
        # Move player
        board[nx][ny] = '*' if board[nx][ny] == '@' else '%'
        # Clear old player position
        board[px][py] = '-' if board[px][py] == '*' else 'X'
        return board
        
    return None

def is_complete(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['X', '%']:
                goals += 1
            if cell in ['$', '%']:
                boxes_on_goals += 1
    return goals == boxes_on_goals and '@' not in ''.join(''.join(row) for row in board)

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_sokoban(initial):
    initial_board = [list(row) for row in initial.strip().split('\n')]
    queue = deque([(initial_board, "")])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_complete(current_board):
            return path
            
        for move in ['U', 'D', 'L', 'R']:
            new_board = make_move(current_board, move)
            if new_board:
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, path + move))
    return None

# Test the puzzle
initial = """++++++
+*$--+
+@@-++
+XX+-+
+++$-+
++++++"""

solution = solve_sokoban(initial)
print(solution)