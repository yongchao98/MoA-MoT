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
    return players, boxes, goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def is_valid_move(board, pos):
    rows, cols = len(board), len(board[0])
    return 0 <= pos[0] < rows and 0 <= pos[1] < cols and board[pos[0]][pos[1]] != '+'

def get_state_key(players, boxes):
    return (tuple(sorted(players)), tuple(sorted(boxes)))

def solve_sokoban(board_str):
    board = parse_board(board_str)
    players, boxes, goals = get_positions(board)
    
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    visited = set()
    queue = deque([(players, set(boxes), '')])
    
    while queue:
        curr_players, curr_boxes, path = queue.popleft()
        state_key = get_state_key(curr_players, tuple(curr_boxes))
        
        if state_key in visited:
            continue
        visited.add(state_key)
        
        if all(box in goals for box in curr_boxes):
            return path
        
        for player in curr_players:
            pi, pj = player
            for move, (di, dj) in moves.items():
                new_pi, new_pj = pi + di, pj + dj
                
                if not is_valid_move(board, (new_pi, new_pj)):
                    continue
                
                if (new_pi, new_pj) in curr_boxes:
                    new_box_i, new_box_j = new_pi + di, new_pj + dj
                    if is_valid_move(board, (new_box_i, new_box_j)) and (new_box_i, new_box_j) not in curr_boxes:
                        new_boxes = set(curr_boxes)
                        new_boxes.remove((new_pi, new_pj))
                        new_boxes.add((new_box_i, new_box_j))
                        
                        new_players = list(curr_players)
                        new_players.remove(player)
                        new_players.append((new_pi, new_pj))
                        
                        queue.append((new_players, new_boxes, path + move))
                else:
                    new_players = list(curr_players)
                    new_players.remove(player)
                    new_players.append((new_pi, new_pj))
                    
                    queue.append((new_players, curr_boxes, path + move))
    
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