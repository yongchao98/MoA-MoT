from collections import deque
import heapq
from copy import deepcopy

def parse_board(board):
    player = None
    boxes = set()
    goals = set()
    walls = set()
    rows = len(board)
    cols = len(board[0])
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['*', '%']:
                player = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '%', '$']:
                goals.add((i, j))
            if board[i][j] == '+':
                walls.add((i, j))
    
    return player, boxes, goals, walls

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_heuristic(boxes, goals):
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

def is_valid_push(player, box, new_box, walls, boxes):
    if new_box in walls or new_box in boxes:
        return False
    return True

def get_moves(player, boxes, walls):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for direction, (dy, dx) in directions.items():
        new_pos = (player[0] + dy, player[1] + dx)
        
        if new_pos in walls:
            continue
            
        if new_pos in boxes:
            new_box_pos = (new_pos[0] + dy, new_pos[1] + dx)
            if is_valid_push(player, new_pos, new_box_pos, walls, boxes - {new_pos}):
                moves.append((direction, new_pos, new_box_pos))
        else:
            moves.append((direction, new_pos, None))
    
    return moves

def solve_sokoban(board):
    player, boxes, goals, walls = parse_board([list(row) for row in board.split('\n')])
    
    initial_state = (player, frozenset(boxes))
    visited = set([initial_state])
    pq = [(get_heuristic(boxes, goals), 0, '', initial_state)]
    
    while pq:
        _, cost, path, (player, boxes) = heapq.heappop(pq)
        
        if all(box in goals for box in boxes):
            return path
        
        for move, new_player, new_box in get_moves(player, boxes, walls):
            new_boxes = set(boxes)
            if new_box:
                new_boxes.remove(new_player)
                new_boxes.add(new_box)
            new_boxes = frozenset(new_boxes)
            
            new_state = (new_player, new_boxes)
            if new_state not in visited:
                visited.add(new_state)
                new_path = path + move
                new_cost = cost + 1
                heapq.heappush(pq, (new_cost + get_heuristic(new_boxes, goals),
                                   new_cost, new_path, new_state))
    return None

# Test the solver
board = """+ + + + + + +
+ - - - @ X +
+ - - - @ - +
+ - @ - X - +
+ - X - @ - +
+ - - % $ + +
+ + + + + + +"""

solution = solve_sokoban(board)
print(solution)