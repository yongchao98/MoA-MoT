from collections import deque
import heapq

def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def find_player(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def find_boxes_and_goals(board):
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return boxes, goals

def is_valid_move(board, pos):
    return 0 <= pos[0] < len(board) and 0 <= pos[1] < len(board[0]) and board[pos[0]][pos[1]] != '+'

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(boxes, goals):
    if not boxes or not goals:
        return 0
    total = 0
    unmatched_boxes = list(boxes)
    unmatched_goals = list(goals)
    while unmatched_boxes:
        min_dist = float('inf')
        min_box = None
        min_goal = None
        for box in unmatched_boxes:
            for goal in unmatched_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    min_box = box
                    min_goal = goal
        total += min_dist
        unmatched_boxes.remove(min_box)
        unmatched_goals.remove(min_goal)
    return total

def get_neighbors(pos):
    return [(pos[0]+1, pos[1]), (pos[0]-1, pos[1]), (pos[0], pos[1]+1), (pos[0], pos[1]-1)]

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player = find_player(board)
    boxes, goals = find_boxes_and_goals(board)
    
    initial_state = (player, frozenset(boxes))
    visited = set([initial_state])
    pq = [(heuristic(boxes, goals), 0, "", initial_state)]
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    while pq:
        _, cost, path, state = heapq.heappop(pq)
        player, boxes = state
        
        if boxes == goals:
            return path
            
        for move, (dy, dx) in moves.items():
            new_player = (player[0] + dy, player[1] + dx)
            
            if not is_valid_move(board, new_player):
                continue
                
            if new_player in boxes:
                new_box_pos = (new_player[0] + dy, new_player[1] + dx)
                if not is_valid_move(board, new_box_pos) or new_box_pos in boxes:
                    continue
                    
                new_boxes = frozenset(b if b != new_player else new_box_pos for b in boxes)
                new_state = (new_player, new_boxes)
                
                if new_state not in visited:
                    visited.add(new_state)
                    new_cost = cost + 1
                    priority = new_cost + heuristic(new_boxes, goals)
                    heapq.heappush(pq, (priority, new_cost, path + move, new_state))
            else:
                new_state = (new_player, boxes)
                if new_state not in visited:
                    visited.add(new_state)
                    new_cost = cost + 1
                    priority = new_cost + heuristic(boxes, goals)
                    heapq.heappush(pq, (priority, new_cost, path + move, new_state))
    
    return None

board = """+ + + + + + + + +  
+ $ + + - X - * +  
+ + + + - - - - +  
+ - - X - @ X - +  
+ - - X - @ @ - +  
+ - - $ - - - - +  
+ - $ $ - @ - - +  
+ - - - - - - - +  
+ + + + + + + + +"""

solution = solve_sokoban(board)
print(solution)