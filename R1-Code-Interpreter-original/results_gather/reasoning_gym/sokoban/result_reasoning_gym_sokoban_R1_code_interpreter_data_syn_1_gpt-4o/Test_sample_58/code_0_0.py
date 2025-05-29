from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 1),
    'boxes': {(2, 2), (2, 3), (4, 2)},
    'goals': {(2, 1), (3, 1), (4, 3)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 6),
              (2, 0), (2, 6),
              (3, 0), (3, 6),
              (4, 0), (4, 6),
              (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6)}
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
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

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state['player'], initial_state['boxes'], "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        player, boxes, path = queue.popleft()
        
        if is_solved(boxes, initial_state['goals']):
            return path
        
        for move, direction in moves.items():
            if is_valid_move(player, direction, boxes, initial_state['walls']):
                new_player, new_boxes = apply_move(player, direction, boxes)
                state = (new_player, frozenset(new_boxes))
                if state not in visited:
                    visited.add(state)
                    queue.append((new_player, new_boxes, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)