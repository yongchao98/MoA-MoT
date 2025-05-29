from collections import deque

def parse_board(board_str):
    return [list(line.strip().replace(" ", "")) for line in board_str.strip().split('\n')]

def find_elements(board):
    player = None
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '%']:
                player = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return player, boxes, goals

def is_valid_move(board, pos):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player, boxes, goals = find_elements(board)
    
    directions = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    
    visited = set()
    queue = deque([(player, frozenset(boxes), "")])
    visited.add((player, frozenset(boxes)))
    
    max_moves = 30  # Limit the search depth
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        if set(curr_boxes) == goals:
            return path
        
        curr_boxes = set(curr_boxes)
        for move, (dy, dx) in directions.items():
            new_player = (curr_player[0] + dy, curr_player[1] + dx)
            
            if not is_valid_move(board, new_player):
                continue
                
            new_boxes = set(curr_boxes)
            if new_player in curr_boxes:
                new_box_pos = (new_player[0] + dy, new_player[1] + dx)
                
                if not is_valid_move(board, new_box_pos) or new_box_pos in curr_boxes:
                    continue
                    
                new_boxes.remove(new_player)
                new_boxes.add(new_box_pos)
            
            new_state = (new_player, frozenset(new_boxes))
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_player, frozenset(new_boxes), path + move))
    
    return None

board = """+++++++++++
+X--+-.-$+
+-@---.-++
+--X--X%++
+-@--$@@@$+
+---------+
+++++++++++"""

solution = solve_sokoban(board)
print(solution)