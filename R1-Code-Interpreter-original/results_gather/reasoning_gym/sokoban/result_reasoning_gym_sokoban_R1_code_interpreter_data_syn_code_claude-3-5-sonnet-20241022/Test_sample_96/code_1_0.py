def parse_board(board_str):
    # Remove spaces between characters
    return [list(row.replace(" ", "")) for row in board_str.strip().split('\n')]

def find_player(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def find_boxes(board):
    boxes = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
    return boxes

def find_goals(board):
    goals = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['X', '$', '%']:
                goals.append((i, j))
    return goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_state_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, player_pos, direction):
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    if not (0 <= new_pos[0] < len(board) and 0 <= new_pos[1] < len(board[0])):
        return False
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if not (0 <= box_pos[0] < len(board) and 0 <= box_pos[1] < len(board[0])):
            return False
        if board[box_pos[0]][box_pos[1]] in ['+', '@', '$']:
            return False
    return True

def make_move(board, player_pos, direction):
    new_board = [row[:] for row in board]
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    # Update player position
    if board[player_pos[0]][player_pos[1]] == '%':
        new_board[player_pos[0]][player_pos[1]] = 'X'
    else:
        new_board[player_pos[0]][player_pos[1]] = '-'
    
    if board[new_pos[0]][new_pos[1]] == 'X':
        new_board[new_pos[0]][new_pos[1]] = '%'
    else:
        new_board[new_pos[0]][new_pos[1]] = '*'
    
    # If pushing a box
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[new_pos[0]][new_pos[1]] == '$':
            new_board[new_pos[0]][new_pos[1]] = '%'
        else:
            new_board[new_pos[0]][new_pos[1]] = '*'
            
        if board[box_pos[0]][box_pos[1]] == 'X':
            new_board[box_pos[0]][box_pos[1]] = '$'
        else:
            new_board[box_pos[0]][box_pos[1]] = '@'
    
    return new_board

def is_solved(board):
    for row in board:
        for cell in row:
            if cell in ['@', 'X']:
                return False
    return True

def solve_sokoban(initial_board):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    start_pos = find_player(initial_board)
    goals = find_goals(initial_board)
    
    queue = [(0, '', initial_board)]
    visited = set()
    
    while queue:
        _, path, current_board = heapq.heappop(queue)
        
        if is_solved(current_board):
            return path
        
        state_string = get_state_string(current_board)
        if state_string in visited:
            continue
        visited.add(state_string)
        
        player_pos = find_player(current_board)
        
        for move, dir_vec in directions.items():
            if is_valid_move(current_board, player_pos, dir_vec):
                new_board = make_move(current_board, player_pos, dir_vec)
                new_boxes = find_boxes(new_board)
                
                # Calculate heuristic
                heuristic = 0
                for box in new_boxes:
                    min_dist = min(manhattan_distance(box, goal) for goal in goals)
                    heuristic += min_dist
                
                heapq.heappush(queue, (len(path) + 1 + heuristic, path + move, new_board))
    
    return None

# Initial board
board_str = """+++++++++
++$%@-+
+++@$-+
+$+---+
++----+
++----+
++--@X+
+X----+
+++++++++"""

import heapq
solution = solve_sokoban(parse_board(board_str))
print(solution)