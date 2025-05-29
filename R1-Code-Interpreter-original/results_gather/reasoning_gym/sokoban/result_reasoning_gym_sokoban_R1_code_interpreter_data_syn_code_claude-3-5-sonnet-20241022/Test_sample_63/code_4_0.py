from collections import deque

def get_player_pos(board):
    for y in range(len(board)):
        for x in range(len(board[y])):
            if board[y][x] in ['*', '%']:
                return (y, x)
    return None

def make_move(board, move):
    board = [list(row) for row in board]
    player = get_player_pos(board)
    if not player:
        return None
    
    py, px = player
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dy, dx = moves[move]
    
    new_py, new_px = py + dy, px + dx
    
    # Check boundaries and walls
    if new_py < 0 or new_px < 0 or new_py >= len(board) or new_px >= len(board[0]) or board[new_py][new_px] == '+':
        return None

    # Moving onto empty space or goal
    if board[new_py][new_px] in ['-', 'X']:
        on_goal = board[py][px] == '%'
        board[py][px] = 'X' if on_goal else '-'
        board[new_py][new_px] = '%' if board[new_py][new_px] == 'X' else '*'
    
    # Moving a box
    elif board[new_py][new_px] in ['@', '$']:
        box_py, box_px = new_py + dy, new_px + dx
        
        # Check if box can be moved
        if (box_py < 0 or box_px < 0 or box_py >= len(board) or box_px >= len(board[0]) or 
            board[box_py][box_px] in ['+', '@', '$']):
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
    boxes = 0
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['@', '$']:
                boxes += 1
            if cell in ['X', '$', '%']:
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return boxes == goals == boxes_on_goals

def solve_sokoban(board):
    queue = deque([(board, '')])
    seen = {tuple(board)}
    max_moves = 50  # Reasonable limit for this puzzle
    
    while queue:
        current, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        if is_solved(current):
            return path
            
        for move in 'UDLR':
            next_state = make_move(current, move)
            if next_state:
                state_tuple = tuple(next_state)
                if state_tuple not in seen:
                    seen.add(state_tuple)
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