from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (3, 4),
    'boxes': {(2, 3), (2, 5), (4, 2), (5, 1), (5, 4)},
    'goals': {(1, 6), (2, 4), (3, 3), (3, 5), (4, 3), (4, 4), (4, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
              (1, 0), (1, 7), (2, 0), (2, 7), (3, 0), (3, 7), (4, 0), (4, 7),
              (5, 0), (5, 7), (6, 0), (6, 7), (7, 0), (7, 1), (7, 2), (7, 3),
              (7, 4), (7, 5), (7, 6), (7, 7)}
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_valid_move(player, direction, boxes, walls):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    if new_player in walls:
        return False
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        if new_box in walls or new_box in boxes:
            return False
    return True

def move(state, move):
    player, boxes = state['player'], state['boxes']
    direction = moves[move]
    new_player = (player[0] + direction[0], player[1] + direction[1])
    new_boxes = set(boxes)
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    return {'player': new_player, 'boxes': new_boxes}

def is_solved(state):
    return state['boxes'] == state['goals']

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state['player'], moves[move], current_state['boxes'], initial_state['walls']):
                new_state = move(current_state, move)
                state_key = (new_state['player'], frozenset(new_state['boxes']))
                if state_key not in visited:
                    visited.add(state_key)
                    queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)