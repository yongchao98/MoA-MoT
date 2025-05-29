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
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + dir[0], new_pos[1] + dir[1])
        if board[box_new_pos[0]][box_new_pos[1]] == '+' or box_new_pos in boxes:
            return False
    return True

def dfs(board, pos, boxes, goals, path, depth, max_depth, visited):
    if depth > max_depth:
        return None
    if boxes == goals:
        return path
    
    state = (pos, frozenset(boxes))
    if state in visited:
        return None
    visited.add(state)
    
    directions = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    
    for move, dir in directions.items():
        if is_valid_move(board, pos, dir, boxes):
            new_pos = (pos[0] + dir[0], pos[1] + dir[1])
            new_boxes = set(boxes)
            
            if new_pos in boxes:
                box_new_pos = (new_pos[0] + dir[0], new_pos[1] + dir[1])
                new_boxes.remove(new_pos)
                new_boxes.add(box_new_pos)
            
            result = dfs(board, new_pos, new_boxes, goals, path + move, depth + 1, max_depth, visited)
            if result:
                return result
    
    return None

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player, boxes, goals = find_elements(board)
    
    # Try increasingly deeper searches
    for max_depth in range(1, 21):  # Limit to depth 20
        visited = set()
        result = dfs(board, player, boxes, goals, "", 0, max_depth, visited)
        if result:
            return result
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