from collections import deque
import heapq
import copy

def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def get_positions(board):
    player = None
    boxes = set()
    goals = set()
    rows, cols = len(board), len(board[0])
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['@', '*']:
                player = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '*']:
                goals.add((i, j))
    
    return player, boxes, goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_heuristic(boxes, goals):
    if not boxes or not goals:
        return 0
    total = 0
    boxes = list(boxes)
    goals = list(goals)
    for box in boxes:
        min_dist = min(manhattan_distance(box, goal) for goal in goals)
        total += min_dist
    return total

def is_valid_move(board, player, boxes, new_player_pos, direction):
    if board[new_player_pos[0]][new_player_pos[1]] == '+':
        return False, None
    
    new_box_pos = None
    if new_player_pos in boxes:
        dx, dy = direction
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if board[new_box_pos[0]][new_box_pos[1]] == '+' or new_box_pos in boxes:
            return False, None
    
    return True, new_box_pos

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player, boxes, goals = get_positions(board)
    
    directions = {
        'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)
    }
    
    initial_state = (player, frozenset(boxes))
    visited = set([initial_state])
    pq = [(get_heuristic(boxes, goals), 0, '', initial_state)]
    
    while pq:
        _, cost, path, (curr_player, curr_boxes) = heapq.heappop(pq)
        
        if curr_boxes == goals:
            return path
        
        for move, (dx, dy) in directions.items():
            new_player_pos = (curr_player[0] + dx, curr_player[1] + dy)
            is_valid, new_box_pos = is_valid_move(board, curr_player, curr_boxes, new_player_pos, (dx, dy))
            
            if not is_valid:
                continue
            
            new_boxes = set(curr_boxes)
            if new_box_pos:
                new_boxes.remove(new_player_pos)
                new_boxes.add(new_box_pos)
            
            new_state = (new_player_pos, frozenset(new_boxes))
            if new_state not in visited:
                visited.add(new_state)
                new_path = path + move
                h_score = get_heuristic(new_boxes, goals)
                heapq.heappush(pq, (h_score + len(new_path), len(new_path), new_path, new_state))
    
    return None

# Test the puzzle
puzzle = """
+ + + + + + + + +
+ - - - + X - $ +
+ X - - + $ - * +
+ - - @ - - X - +
+ @ $ - - - @ - +
+ - - - - - - - +
+ + + + + + + + +
"""

solution = solve_sokoban(puzzle)
print(solution)