def parse_board(board_str):
    # Remove any whitespace and split into lines
    lines = [line.strip().replace(" ", "") for line in board_str.strip().split('\n')]
    return [list(line) for line in lines]

def find_player(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '%']:
                return (i, j)
    return None

def find_boxes(board):
    boxes = set()
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def find_goals(board):
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return goals

def is_valid_move(board, player, direction):
    new_pos = (player[0] + direction[0], player[1] + direction[1])
    if (new_pos[0] < 0 or new_pos[0] >= len(board) or 
        new_pos[1] < 0 or new_pos[1] >= len(board[0])):
        return False
    return board[new_pos[0]][new_pos[1]] != '+'

def make_move(board, player, boxes, direction):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    new_boxes = set(boxes)
    
    # Check if pushing a box
    if new_player in boxes:
        new_box_pos = (new_player[0] + direction[0], new_player[1] + direction[1])
        # Check if new box position is valid
        if (new_box_pos[0] < 0 or new_box_pos[0] >= len(board) or 
            new_box_pos[1] < 0 or new_box_pos[1] >= len(board[0]) or
            new_box_pos in boxes or board[new_box_pos[0]][new_box_pos[1]] == '+'):
            return None, None
        new_boxes.remove(new_player)
        new_boxes.add(new_box_pos)
    
    return new_player, new_boxes

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(boxes, goals):
    if not boxes or not goals:
        return 0
    total = 0
    for box in boxes:
        min_dist = min(manhattan_distance(box, goal) for goal in goals)
        total += min_dist
    return total

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player = find_player(board)
    boxes = find_boxes(board)
    goals = find_goals(board)
    
    if not player or not boxes or not goals:
        return None
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    import heapq
    start_state = (player, frozenset(boxes))
    visited = set([start_state])
    pq = [(heuristic(boxes, goals), 0, '', start_state)]
    
    while pq:
        _, cost, path, (curr_player, curr_boxes) = heapq.heappop(pq)
        
        if curr_boxes == goals:
            return path
        
        for move, dir_vec in directions.items():
            if not is_valid_move(board, curr_player, dir_vec):
                continue
                
            new_player, new_boxes = make_move(board, curr_player, curr_boxes, dir_vec)
            if new_player is None:
                continue
                
            new_state = (new_player, frozenset(new_boxes))
            if new_state in visited:
                continue
                
            visited.add(new_state)
            new_cost = cost + 1
            priority = new_cost + heuristic(new_boxes, goals)
            heapq.heappush(pq, (priority, new_cost, path + move, new_state))
    
    return None

# Input board - exactly as provided
board = """+++++++++++
+X--+-.-$+
+-@---.-++
+--X--X%++
+-@--$@@@$+
+---------+
+++++++++++"""

solution = solve_sokoban(board)
print(solution)