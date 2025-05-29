from heapq import heappush, heappop
import numpy as np

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def calculate_heuristic(boxes, goals):
    if not boxes:
        return 0
    total = 0
    used_goals = set()
    for box in boxes:
        min_dist = float('inf')
        best_goal = None
        for goal in goals:
            if goal not in used_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    best_goal = goal
        if best_goal:
            used_goals.add(best_goal)
            total += min_dist
    return total

def get_next_moves(board, player_pos, boxes):
    moves = []
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    
    for direction, dy, dx in directions:
        new_y, new_x = player_pos[0] + dy, player_pos[1] + dx
        
        if board[new_y][new_x] == '+':
            continue
            
        if (new_y, new_x) in boxes:
            box_new_y, box_new_x = new_y + dy, new_x + dx
            if board[box_new_y][box_new_x] != '+' and (box_new_y, box_new_x) not in boxes:
                new_boxes = set(boxes)
                new_boxes.remove((new_y, new_x))
                new_boxes.add((box_new_y, box_new_x))
                moves.append((direction, (new_y, new_x), frozenset(new_boxes)))
        else:
            moves.append((direction, (new_y, new_x), frozenset(boxes)))
    
    return moves

def solve_sokoban(board):
    # Find initial state
    boxes = set()
    goals = set()
    player_pos = None
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
            if board[i][j] in ['*', '%']:
                player_pos = (i, j)
    
    # A* search
    start_state = (player_pos, frozenset(boxes))
    visited = {start_state}
    pq = [(calculate_heuristic(boxes, goals), 0, '', start_state)]
    max_steps = 100  # Limit search depth
    
    while pq:
        _, cost, path, (curr_pos, curr_boxes) = heappop(pq)
        
        if len(path) > max_steps:
            continue
            
        if set(curr_boxes) == goals:
            return path
            
        for direction, new_pos, new_boxes in get_next_moves(board, curr_pos, curr_boxes):
            new_state = (new_pos, new_boxes)
            if new_state not in visited:
                visited.add(new_state)
                new_cost = cost + 1
                h_cost = calculate_heuristic(new_boxes, goals)
                heappush(pq, (new_cost + h_cost, new_cost, path + direction, new_state))
    
    return None

# Parse input
board = []
puzzle = [
    "+ + + + + + + +",
    "+ * @ - - - X +",
    "+ @ $ + - @ - +",
    "+ X + X - - - +",
    "+ + + X - X - +",
    "+ + X - - - - +",
    "+ - - @ X @ - +",
    "+ - @ - - @ - +",
    "+ - - - - - - +",
    "+ + + + + + + +"
]
board = [''.join(row.split()) for row in puzzle]

solution = solve_sokoban(board)
print(solution)