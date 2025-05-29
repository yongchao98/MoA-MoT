from collections import deque

def get_state(board):
    # Convert board to a string for state comparison
    return ''.join(''.join(row) for row in board)

def find_player(board):
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell in ['*', '%']:
                return i, j
    return None

def try_move(board, move):
    board = [list(row) for row in board]
    px, py = find_player(board)
    
    # Define move directions
    if move == 'U': dx, dy = -1, 0
    elif move == 'D': dx, dy = 1, 0
    elif move == 'L': dx, dy = 0, -1
    elif move == 'R': dx, dy = 0, 1
    
    # New position
    nx, ny = px + dx, py + dy
    
    # Check bounds and walls
    if nx < 0 or ny < 0 or nx >= len(board) or ny >= len(board[0]) or board[nx][ny] == '+':
        return None
    
    # Moving to empty space or goal
    if board[nx][ny] in ['-', 'X']:
        # Update player position
        new_player = '*' if board[nx][ny] == '-' else '%'
        old_space = '-' if board[px][py] == '*' else 'X'
        board[nx][ny] = new_player
        board[px][py] = old_space
        return board
    
    # Moving box
    if board[nx][ny] in ['@', '$']:
        # Check box destination
        bx, by = nx + dx, ny + dy
        if (bx < 0 or by < 0 or bx >= len(board) or by >= len(board[0]) or 
            board[bx][by] in ['+', '@', '$']):
            return None
        
        # Move box
        new_box = '@' if board[bx][by] == '-' else '$'
        board[bx][by] = new_box
        
        # Move player
        new_player = '*' if board[nx][ny] == '@' else '%'
        board[nx][ny] = new_player
        
        # Clear old player position
        board[px][py] = '-' if board[px][py] == '*' else 'X'
        
        return board
    
    return None

def is_goal(board):
    # Check if all boxes are on goals
    boxes = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['@', '$']:
                boxes += 1
            if cell == '$':
                boxes_on_goals += 1
    return boxes == boxes_on_goals

def solve(initial):
    initial_board = [list(row) for row in initial.strip().split('\n')]
    queue = deque([(initial_board, "")])
    seen = {get_state(initial_board)}
    
    while queue:
        board, path = queue.popleft()
        
        if is_goal(board):
            return path
        
        for move in ['U', 'D', 'L', 'R']:
            new_board = try_move(board, move)
            if new_board:
                state = get_state(new_board)
                if state not in seen:
                    seen.add(state)
                    queue.append((new_board, path + move))
    return "No solution found"

# Initial puzzle
puzzle = """++++++
+*$--+
+@@-++
+XX+-+
+++$-+
++++++"""

result = solve(puzzle)
print(result)