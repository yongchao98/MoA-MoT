from collections import deque

def get_player_pos(board):
    for y in range(len(board)):
        for x in range(len(board[y])):
            if board[y][x] in ['*', '%']:
                return (y, x)
    return None

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_goals_and_boxes(board):
    goals = []
    boxes = []
    for y in range(len(board)):
        for x in range(len(board[y])):
            if board[y][x] in ['X', '$', '%']:
                goals.append((y, x))
            if board[y][x] in ['@', '$']:
                boxes.append((y, x))
    return goals, boxes

def heuristic(board):
    goals, boxes = get_goals_and_boxes(board)
    if len(boxes) != len(goals):
        return float('inf')
    
    total_dist = 0
    for box in boxes:
        min_dist = min(manhattan_distance(box, goal) for goal in goals)
        total_dist += min_dist
    return total_dist

def make_move(board, move):
    board = [list(row) for row in board]
    player = get_player_pos(board)
    if not player:
        return None
    
    py, px = player
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dy, dx = moves[move]
    
    new_py, new_px = py + dy, px + dx
    
    if board[new_py][new_px] == '+':
        return None
        
    if board[new_py][new_px] in ['-', 'X']:
        on_goal = board[py][px] == '%'
        board[py][px] = 'X' if on_goal else '-'
        board[new_py][new_px] = '%' if board[new_py][new_px] == 'X' else '*'
    
    elif board[new_py][new_px] in ['@', '$']:
        box_py, box_px = new_py + dy, new_px + dx
        
        if board[box_py][box_px] in ['+', '@', '$']:
            return None
            
        on_goal = board[new_py][new_px] == '$'
        new_box_on_goal = board[box_py][box_px] == 'X'
        board[box_py][box_px] = '$' if new_box_on_goal else '@'
        
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
    queue = [(heuristic(board), 0, initial, '')]  # (priority, moves, state, path)
    seen = {initial}
    max_moves = 30  # Limit search depth
    
    while queue:
        queue.sort()  # Sort by priority
        _, moves, current, path = queue.pop(0)
        
        if is_solved(current):
            return path
            
        if moves >= max_moves:
            continue
            
        for move in 'UDLR':
            next_state = make_move(current, move)
            if next_state:
                next_state_tuple = tuple(next_state)
                if next_state_tuple not in seen:
                    seen.add(next_state_tuple)
                    priority = heuristic(next_state) + moves + 1
                    queue.append((priority, moves + 1, next_state_tuple, path + move))
                    
                    if len(queue) > 1000:  # Limit queue size
                        queue.sort()
                        queue = queue[:500]
    
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