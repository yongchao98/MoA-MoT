def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '*']:
                return (i, j)
    return None

def get_boxes_and_goals(board):
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['$', '@']:
                boxes.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
            if board[i][j] in ['X', '*']:
                goals.add((i, j))
    return boxes, goals

def is_valid_pos(pos, board):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def get_state_string(board):
    return ''.join(''.join(row) for row in board)

def make_move(board, player_pos, direction, boxes):
    new_board = [row[:] for row in board]
    new_boxes = set(boxes)
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    # Moving box
    if (new_pos[0], new_pos[1]) in boxes:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        new_boxes.remove(new_pos)
        new_boxes.add(box_new_pos)
        
        # Update box position
        if board[box_new_pos[0]][box_new_pos[1]] == 'X':
            new_board[box_new_pos[0]][box_new_pos[1]] = '$'
        else:
            new_board[box_new_pos[0]][box_new_pos[1]] = '@'
    
    # Update old player position
    if board[player_pos[0]][player_pos[1]] == '*':
        new_board[player_pos[0]][player_pos[1]] = 'X'
    else:
        new_board[player_pos[0]][player_pos[1]] = '-'
    
    # Update new player position
    if board[new_pos[0]][new_pos[1]] == 'X':
        new_board[new_pos[0]][new_pos[1]] = '*'
    else:
        new_board[new_pos[0]][new_pos[1]] = '@'
    
    return new_board, new_boxes

def is_valid_move(board, player_pos, direction, boxes):
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    if not is_valid_pos(new_pos, board):
        return False
        
    # If moving to a box, check if box can be pushed
    if (new_pos[0], new_pos[1]) in boxes:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        return (is_valid_pos(box_new_pos, board) and 
                box_new_pos not in boxes)
    
    return True

def solve_sokoban(board):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    boxes, goals = get_boxes_and_goals(board)
    visited = set()
    
    from collections import deque
    queue = deque([(board, '', get_player_pos(board), boxes)])
    
    while queue:
        current_board, path, player_pos, current_boxes = queue.popleft()
        state = (get_state_string(current_board), player_pos)
        
        if state in visited:
            continue
        visited.add(state)
        
        # Check if won
        if all(box in goals for box in current_boxes):
            return path
        
        # Try all directions
        for move, dir_offset in directions.items():
            if is_valid_move(current_board, player_pos, dir_offset, current_boxes):
                new_board, new_boxes = make_move(current_board, player_pos, dir_offset, current_boxes)
                new_player_pos = (player_pos[0] + dir_offset[0], player_pos[1] + dir_offset[1])
                queue.append((new_board, path + move, new_player_pos, new_boxes))
    
    return None

# Initial board
board_str = """
+ + + + + + + +
+ - - - - - - +
+ - $ @ - $ - +
+ - X - @ $ - +
+ - - - X * - +
+ + + + + + + +
"""

initial_board = parse_board(board_str)
solution = solve_sokoban(initial_board)
print(solution)