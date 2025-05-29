from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (4, 4),
    'boxes': {(2, 1), (2, 3), (3, 3), (5, 3), (6, 1)},
    'goals': {(3, 1), (4, 3), (5, 2), (6, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 5),
              (2, 0), (2, 5),
              (3, 0), (3, 5),
              (4, 0), (4, 5),
              (5, 0), (5, 5),
              (6, 0), (6, 5),
              (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5)}
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Check if the move is valid
def is_valid_move(player, direction, boxes, walls):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    if new_player in walls:
        return False
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        if new_box in walls or new_box in boxes:
            return False
    return True

# Apply the move
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

# Check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Perform BFS to find the solution
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
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)