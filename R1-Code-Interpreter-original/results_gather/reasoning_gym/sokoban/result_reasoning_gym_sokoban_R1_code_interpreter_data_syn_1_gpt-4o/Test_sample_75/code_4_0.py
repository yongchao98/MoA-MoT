from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (7, 5),
    'boxes': {(6, 4), (8, 4), (8, 6)},
    'goals': {(5, 5), (6, 3), (8, 3)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
              (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6),
              (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6),
              (4, 0), (4, 6),
              (5, 0), (5, 6),
              (6, 0), (6, 6),
              (7, 0), (7, 6),
              (8, 0), (8, 6),
              (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6)},
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
def apply_move(state, move):
    player = state['player']
    boxes = set(state['boxes'])
    direction = moves[move]
    new_player = (player[0] + direction[0], player[1] + direction[1])
    if new_player in boxes:
        boxes.remove(new_player)
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        boxes.add(new_box)
    return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state['player'], moves[move], current_state['boxes'], current_state['walls']):
                new_state = apply_move(current_state, move)
                state_key = (new_state['player'], frozenset(new_state['boxes']))
                if state_key not in visited:
                    visited.add(state_key)
                    queue.append((new_state, path + move))
                    print(f"Move: {move}, Path: {path + move}, Player: {new_state['player']}, Boxes: {new_state['boxes']}")
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)