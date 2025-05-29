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

def is_deadlock(board, box_pos, goals):
    # Check if box is in corner
    if (board[box_pos[0]-1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]-1] == '+') or \
       (board[box_pos[0]-1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]+1] == '+') or \
       (board[box_pos[0]+1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]-1] == '+') or \
       (board[box_pos[0]+1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]+1] == '+'):
        return box_pos not in goals
    return False

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def min_matching_cost(boxes, goals):
    if len(boxes) != len(goals):
        return float('inf')
    
    unmatched_boxes = list(boxes)
    unmatched_goals = list(goals)
    total_cost = 0
    
    while unmatched_boxes:
        min_dist = float('inf')
        best_pair = None
        for i, box in enumerate(unmatched_boxes):
            for j, goal in enumerate(unmatched_goals):
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    best_pair = (i, j)
        
        if best_pair is None:
            break
            
        box_idx, goal_idx = best_pair
        total_cost += manhattan_distance(unmatched_boxes[box_idx], unmatched_goals[goal_idx])
        unmatched_boxes.pop(box_idx)
        unmatched_goals.pop(goal_idx)
    
    return total_cost

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player, boxes, goals = find_elements(board)
    
    if not player or not boxes or not goals:
        return None
    
    directions = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    
    import heapq
    from collections import defaultdict
    
    visited = set()
    min_cost = defaultdict(lambda: float('inf'))
    start_state = (player, frozenset(boxes))
    
    pq = [(min_matching_cost(boxes, goals), 0, '', start_state)]
    
    while pq:
        h_cost, g_cost, path, (curr_player, curr_boxes) = heapq.heappop(pq)
        
        if len(path) > 100:  # Limit solution length
            continue
            
        state_hash = (curr_player, curr_boxes)
        if state_hash in visited and g_cost >= min_cost[state_hash]:
            continue
            
        visited.add(state_hash)
        min_cost[state_hash] = g_cost
        
        if curr_boxes == goals:
            return path
        
        curr_boxes = set(curr_boxes)
        for move, (dy, dx) in directions.items():
            new_player = (curr_player[0] + dy, curr_player[1] + dx)
            
            # Check if move is valid
            if board[new_player[0]][new_player[1]] == '+':
                continue
                
            new_boxes = set(curr_boxes)
            if new_player in curr_boxes:
                new_box_pos = (new_player[0] + dy, new_player[1] + dx)
                
                if (board[new_box_pos[0]][new_box_pos[1]] == '+' or
                    new_box_pos in curr_boxes or
                    is_deadlock(board, new_box_pos, goals)):
                    continue
                    
                new_boxes.remove(new_player)
                new_boxes.add(new_box_pos)
            
            new_state = (new_player, frozenset(new_boxes))
            if new_state not in visited or g_cost + 1 < min_cost[new_state]:
                new_h_cost = min_matching_cost(new_boxes, goals)
                heapq.heappush(pq, (new_h_cost + g_cost + 1, g_cost + 1, path + move, new_state))
    
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