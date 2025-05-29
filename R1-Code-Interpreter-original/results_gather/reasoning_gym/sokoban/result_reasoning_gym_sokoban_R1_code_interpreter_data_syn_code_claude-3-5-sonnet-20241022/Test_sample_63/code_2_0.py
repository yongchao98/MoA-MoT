from collections import deque

def get_player_pos(board):
    for y in range(len(board)):
        for x in range(len(board[y])):
            if board[y][x] in ['*', '%']:
                return (y, x)
    return None

def make_move(board, move):
    # Convert string board to list of lists
    board = [list(row) for row in board]
    player = get_player_pos(board)
    if not player:
        return None
    
    py, px = player
    
    # Define move directions
    if move == 'U': dy, dx = -1, 0
    elif move == 'D': dy, dx = 1, 0
    elif move == 'L': dy, dx = 0, -1
    elif move == 'R': dy, dx = 0, 1
    else: return None
    
    # New position
    new_py, new_px = py + dy, px + dx
    
    # Check if move is valid
    if board[new_py][new_px] == '+':
        return None
        
    # Moving onto empty space or goal
    if board[new_py][new_px] in ['-', 'X']:
        # Update player position
        on_goal = board[py][px] == '%'
        board[py][px] = 'X' if on_goal else '-'
        board[new_py][new_px] = '%' if board[new_py][new_px] == 'X' else '*'
        
    # Moving a box
    elif board[new_py][new_px] in ['@', '$']:
        box_py, box_px = new_py + dy, new_px + dx
        
        # Check if box can be moved
        if board[box_py][box_px] in ['+', '@', '$']:
            return None
            
        # Move box
        on_goal = board[new_py][new_px] == '$'
        new_box_on_goal = board[box_py][box_px] == 'X'
        board[box_py][box_px] = '$' if new_box_on_goal else '@'
        
        # Move player
        player_on_goal = board[py][px] == '%'
        board[py][px] = 'X' if player_on_goal else '-'
        board[new_py][new_px] = '%' if on_goal else '*'
    else:
        return None
        
    return [''.join(row) for row in board]

def is_solved(board):
    for row in board:
        if '@' in row or 'X' in row:
            return False
    return True

def solve_sokoban(board):
    initial = tuple(board)
    queue = deque([(initial, '')])
    seen = {initial}
    
    while queue:
        current, path = queue.popleft()
        
        if is_solved(current):
            return path
            
        for move in 'UDLR':
            next_state = make_move(current, move)
            if next_state:
                next_state = tuple(next_state)
                if next_state not in seen:
                    seen.add(next_state)
                    queue.append((next_state, path + move))
    return None

# Initial board
board = [
    "+++++++++++",
    "+-$$X--@-+",
    "+-----X@-+",
    "+-@X@@-X-+",
    "+-@-----+",
    "+-@-X--$-+",
    "+-*X----+",
    "+-----X--+",
    "+---------+",
    "+++++++++++"
]

solution = solve_sokoban(board)
print(solution if solution else "No solution found")