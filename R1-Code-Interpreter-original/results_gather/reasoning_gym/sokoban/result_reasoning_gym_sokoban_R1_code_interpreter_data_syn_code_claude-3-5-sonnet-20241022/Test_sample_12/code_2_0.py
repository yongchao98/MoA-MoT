def parse_board(board_str):
    return [list(row.strip().replace(" ", "")) for row in board_str.strip().split('\n') if row.strip()]

def find_elements(board):
    player = None
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                player = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return player, boxes, goals

def is_valid_move(board, pos, dir, boxes):
    new_pos = (pos[0] + dir[0], pos[1] + dir[1])
    if not (0 <= new_pos[0] < len(board) and 0 <= new_pos[1] < len(board[0])):
        return False
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + dir[0], new_pos[1] + dir[1])
        if not (0 <= box_new_pos[0] < len(board) and 0 <= box_new_pos[1] < len(board[0])):
            return False
        if board[box_new_pos[0]][box_new_pos[1]] == '+' or box_new_pos in boxes:
            return False
    return True

def simple_heuristic(boxes, goals):
    return sum(min(abs(b[0] - g[0]) + abs(b[1] - g[1]) for g in goals) for b in boxes)

def solve_sokoban(board_str, max_depth=50):
    board = parse_board(board_str)
    player, boxes, goals = find_elements(board)
    
    directions = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    
    from collections import deque
    queue = deque([(player, frozenset(boxes), "", 0)])
    visited = {(player, frozenset(boxes))}
    
    while queue:
        pos, boxes, path, depth = queue.popleft()
        
        if depth >= max_depth:
            continue
            
        if boxes == goals:
            return path
            
        for move, dir in directions.items():
            if is_valid_move(board, pos, dir, boxes):
                new_pos = (pos[0] + dir[0], pos[1] + dir[1])
                new_boxes = set(boxes)
                
                if new_pos in boxes:
                    box_new_pos = (new_pos[0] + dir[0], new_pos[1] + dir[1])
                    new_boxes.remove(new_pos)
                    new_boxes.add(box_new_pos)
                
                new_boxes = frozenset(new_boxes)
                new_state = (new_pos, new_boxes)
                
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_pos, new_boxes, path + move, depth + 1))
    
    return None

puzzle = """
++++++++++ 
+$+++X-@-+
+++++X@@-+
+X++X--X-+
+-+X@--X-+
+@--X@---+
+---@-@--+
+-----*--+
++++++++++ 
"""

solution = solve_sokoban(puzzle)
print(solution)