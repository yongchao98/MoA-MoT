from collections import deque
import heapq

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 1),
    'boxes': {(1, 2), (2, 1), (2, 3), (3, 5), (3, 6), (3, 7), (6, 2), (7, 3)},
    'goals': {(1, 6), (2, 5), (6, 1), (7, 4), (8, 2)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
              (1, 0), (1, 7), (1, 8), (2, 0), (2, 6), (2, 7), (2, 8), (3, 0), (3, 8),
              (4, 0), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8),
              (5, 0), (5, 2), (5, 8), (6, 0), (6, 5), (6, 6), (6, 7), (6, 8), (7, 0),
              (7, 5), (7, 6), (7, 7), (7, 8), (8, 0), (8, 1), (8, 3), (8, 4), (8, 5),
              (8, 6), (8, 7), (8, 8)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Function to check if a move is valid
def is_valid_move(player, direction, boxes, walls):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    if new_player in walls:
        return False
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        if new_box in walls or new_box in boxes:
            return False
    return True

# Function to apply a move
def apply_move(player, direction, boxes):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    new_boxes = set(boxes)
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    return new_player, new_boxes

# Function to check if the puzzle is solved
def is_solved(boxes, goals):
    return boxes == goals

# Heuristic function: sum of Manhattan distances from boxes to nearest goals
def heuristic(boxes, goals):
    return sum(min(abs(bx - gx) + abs(by - gy) for gx, gy in goals) for bx, by in boxes)

# A* search to find the solution
def solve_sokoban(initial_state):
    queue = []
    heapq.heappush(queue, (0, initial_state['player'], initial_state['boxes'], ""))
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        _, player, boxes, path = heapq.heappop(queue)
        
        if is_solved(boxes, initial_state['goals']):
            return path
        
        for move, direction in moves.items():
            if is_valid_move(player, direction, boxes, initial_state['walls']):
                new_player, new_boxes = apply_move(player, direction, boxes)
                state = (new_player, frozenset(new_boxes))
                if state not in visited:
                    visited.add(state)
                    new_cost = len(path) + 1 + heuristic(new_boxes, initial_state['goals'])
                    heapq.heappush(queue, (new_cost, new_player, new_boxes, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)