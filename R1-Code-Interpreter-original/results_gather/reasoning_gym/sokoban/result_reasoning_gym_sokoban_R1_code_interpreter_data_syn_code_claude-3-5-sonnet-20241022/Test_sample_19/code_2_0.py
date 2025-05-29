from collections import deque

def parse_board(board_str):
    # Convert string to 2D list, removing spaces between characters
    return [list(''.join(line.split())) for line in board_str.strip().split('\n')]

def get_player_boxes_goals(board):
    player = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in '@%':
                player = (i, j)
                if board[i][j] == '%':
                    goals.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    return player, boxes, goals

def is_valid(pos, board):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def get_moves(pos, boxes, board):
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    valid_moves = []
    
    for dy, dx, move in directions:
        new_pos = (pos[0] + dy, pos[1] + dx)
        if not is_valid(new_pos, board):
            continue
            
        if new_pos in boxes:
            box_new_pos = (new_pos[0] + dy, new_pos[1] + dx)
            if not is_valid(box_new_pos, board) or box_new_pos in boxes:
                continue
            new_boxes = set(boxes)
            new_boxes.remove(new_pos)
            new_boxes.add(box_new_pos)
            valid_moves.append((new_pos, new_boxes, move))
        else:
            valid_moves.append((new_pos, boxes, move))
    
    return valid_moves

def solve(board_str):
    board = parse_board(board_str)
    player, boxes, goals = get_player_boxes_goals(board)
    
    queue = deque([(player, boxes, "")])
    seen = set()
    
    while queue:
        pos, current_boxes, path = queue.popleft()
        state = (pos, frozenset(current_boxes))
        
        if current_boxes == goals:
            return path
            
        if state in seen:
            continue
        seen.add(state)
        
        for new_pos, new_boxes, move in get_moves(pos, current_boxes, board):
            queue.append((new_pos, new_boxes, path + move))
    
    return None

puzzle = """+ + + + + + +
+ $ + - - X +
+ - - - - - +
+ - - - @ - +
+ - X - @ - +
+ - - - X - +
+ - - @ X X +
+ - - - @ @ +
+ - - - @ % +
+ + + + + + +"""

solution = solve(puzzle)
print(solution)