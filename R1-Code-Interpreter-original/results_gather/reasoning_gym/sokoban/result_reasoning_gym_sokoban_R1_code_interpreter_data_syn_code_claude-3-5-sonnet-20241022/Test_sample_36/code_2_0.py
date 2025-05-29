from collections import deque

def get_initial_state(board):
    player = None
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '*', '%']:  # player
                player = (i, j)
            if board[i][j] in ['@', '$']:  # boxes
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '*', '%']:  # goals
                goals.add((i, j))
    return player, boxes, goals

def is_valid_pos(pos, board):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def get_moves(player, boxes, board, goals):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    valid_moves = []
    
    for move, (dx, dy) in directions.items():
        new_player = (player[0] + dx, player[1] + dy)
        
        if not is_valid_pos(new_player, board):
            continue
            
        if new_player in boxes:  # pushing a box
            new_box = (new_player[0] + dx, new_player[1] + dy)
            if is_valid_pos(new_box, board) and new_box not in boxes:
                new_boxes = {new_box if b == new_player else b for b in boxes}
                valid_moves.append((new_player, new_boxes, move))
        else:  # just moving
            valid_moves.append((new_player, boxes, move))
    
    return valid_moves

def encode_state(player, boxes):
    boxes_tuple = tuple(sorted(boxes))
    return (player, boxes_tuple)

def is_solved(boxes, goals):
    return boxes == goals

def solve_sokoban(board):
    player, boxes, goals = get_initial_state(board)
    initial_state = encode_state(player, boxes)
    
    queue = deque([(player, boxes, "")])
    visited = {initial_state}
    
    max_iterations = 100000  # prevent infinite loops
    iterations = 0
    
    while queue and iterations < max_iterations:
        iterations += 1
        current_player, current_boxes, path = queue.popleft()
        
        if is_solved(current_boxes, goals):
            return path
        
        for new_player, new_boxes, move in get_moves(current_player, current_boxes, board, goals):
            state = encode_state(new_player, new_boxes)
            if state not in visited:
                visited.add(state)
                queue.append((new_player, new_boxes, path + move))
    
    return None

def print_board_state(board, player, boxes, goals):
    board_state = [list(row) for row in board]
    for i, j in goals:
        board_state[i][j] = 'X'
    for i, j in boxes:
        if (i, j) in goals:
            board_state[i][j] = '$'
        else:
            board_state[i][j] = '@'
    if player in goals:
        board_state[player[0]][player[1]] = '*'
    else:
        board_state[player[0]][player[1]] = '@'
    for row in board_state:
        print(''.join(row))

# Initial board
board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '+', '-', '$', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '-', '@', '+', '-', '+'],
    ['+', '$', '$', 'X', '-', '@', '+'],
    ['+', '+', '+', '+', '$', '*', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)