from collections import deque

def parse_board(board_str):
    return [list(row.strip().split()) for row in board_str.strip().split('\n')]

def get_positions(board):
    players = []
    boxes = []
    goals = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                players.append((i, j))
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
            if board[i][j] in ['X', '%', '$']:
                goals.append((i, j))
    return players[0], set(boxes), set(goals)  # Only take first player

def is_valid_pos(board, pos):
    return 0 <= pos[0] < len(board) and 0 <= pos[1] < len(board[0]) and board[pos[0]][pos[1]] != '+'

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player, boxes, goals = get_positions(board)
    
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    visited = set()
    queue = deque([(player, boxes, '')])
    max_depth = 30  # Limit search depth
    
    while queue:
        player, curr_boxes, path = queue.popleft()
        
        if len(path) > max_depth:  # Depth limit check
            continue
            
        state = (player, tuple(sorted(curr_boxes)))
        if state in visited:
            continue
        visited.add(state)
        
        if curr_boxes == goals:
            return path
        
        pi, pj = player
        for move, (di, dj) in moves.items():
            new_pi, new_pj = pi + di, pj + dj
            
            if not is_valid_pos(board, (new_pi, new_pj)):
                continue
            
            if (new_pi, new_pj) in curr_boxes:
                new_box_pos = (new_pi + di, new_pj + dj)
                if is_valid_pos(board, new_box_pos) and new_box_pos not in curr_boxes:
                    new_boxes = set(curr_boxes)
                    new_boxes.remove((new_pi, new_pj))
                    new_boxes.add(new_box_pos)
                    queue.append(((new_pi, new_pj), new_boxes, path + move))
            else:
                queue.append(((new_pi, new_pj), curr_boxes, path + move))
    
    return None

board = """+ + + + + + + + + +
+ - - - - - - - - +
+ - @ $ - $ - @ - +
+ - * - - X - - - +
+ X - - - - - @ - +
+ - $ X - - - $ - +
+ - - X @ - X @ - +
+ - + + - - - - - +
+ + + + + + + + + +"""

solution = solve_sokoban(board)
print(solution)