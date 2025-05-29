from heapq import heappush, heappop
import numpy as np

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def parse_board(board_str):
    board = [list(row.strip().split()) for row in board_str.strip().split('\n')]
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == '*':
                player_pos = (i, j)
            elif board[i][j] == '@':
                boxes.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
    
    return board, player_pos, boxes, goals

def is_deadlock(board, box_pos, goals):
    # Check if box is in corner
    if (board[box_pos[0]-1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]-1] == '+') or \
       (board[box_pos[0]-1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]+1] == '+') or \
       (board[box_pos[0]+1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]-1] == '+') or \
       (board[box_pos[0]+1][box_pos[1]] == '+' and board[box_pos[0]][box_pos[1]+1] == '+'):
        return box_pos not in goals
    return False

def heuristic(boxes, goals):
    # Minimum sum of manhattan distances from each box to a goal
    total = 0
    used_goals = set()
    sorted_boxes = sorted(boxes, key=lambda box: min(manhattan_distance(box, goal) for goal in goals))
    
    for box in sorted_boxes:
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

def solve_sokoban(board_str):
    board, initial_player, initial_boxes, goals = parse_board(board_str)
    visited = set()
    pq = [(0, 0, initial_player, initial_boxes, "")]  # (priority, moves, player_pos, boxes, path)
    moves_count = 0
    
    while pq and moves_count < 1000:  # Add move limit to prevent infinite loops
        _, moves, player_pos, boxes, path = heappop(pq)
        state = (player_pos, frozenset(boxes))
        
        if boxes == goals:
            return path
            
        if state in visited:
            continue
            
        visited.add(state)
        moves_count += 1
        
        for direction, (dx, dy) in {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}.items():
            new_pos = (player_pos[0] + dx, player_pos[1] + dy)
            
            # Check boundaries and walls
            if (new_pos[0] < 0 or new_pos[0] >= len(board) or 
                new_pos[1] < 0 or new_pos[1] >= len(board[0]) or
                board[new_pos[0]][new_pos[1]] == '+'):
                continue
                
            new_boxes = set(boxes)
            if new_pos in boxes:
                box_new_pos = (new_pos[0] + dx, new_pos[1] + dy)
                
                # Check if box move is valid
                if (box_new_pos[0] < 0 or box_new_pos[0] >= len(board) or 
                    box_new_pos[1] < 0 or box_new_pos[1] >= len(board[0]) or
                    board[box_new_pos[0]][box_new_pos[1]] == '+' or
                    box_new_pos in boxes):
                    continue
                    
                # Check for deadlock
                if is_deadlock(board, box_new_pos, goals):
                    continue
                    
                new_boxes.remove(new_pos)
                new_boxes.add(box_new_pos)
            
            # Calculate priority using heuristic
            h = heuristic(new_boxes, goals)
            priority = moves + 1 + h
            
            heappush(pq, (priority, moves + 1, new_pos, new_boxes, path + direction))
    
    return None

puzzle = """
+ + + + + + + +
+ + - X - - - +
+ X - @ X @ * +
+ X @ - - - - +
+ + - - - @ - +
+ - - - - @ - +
+ X @ X - - - +
+ + + + + + + +
"""

solution = solve_sokoban(puzzle)
print(f"<<<{solution}>>>")