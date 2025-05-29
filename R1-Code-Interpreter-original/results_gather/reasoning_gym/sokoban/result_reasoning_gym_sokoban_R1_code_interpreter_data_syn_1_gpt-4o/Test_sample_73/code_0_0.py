from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (6, 8),
    'boxes': {(2, 3), (2, 6), (5, 7), (5, 8)},
    'goals': {(3, 4), (4, 8), (6, 4), (6, 7)},
    'walls': {(0, i) for i in range(10)} | {(7, i) for i in range(10)} |
             {(i, 0) for i in range(8)} | {(i, 9) for i in range(8)} |
             {(2, 7), (2, 8), (3, 3), (3, 5), (3, 6), (3, 7), (3, 8), (4, 3), (4, 5), (4, 6), (4, 7), (5, 3), (5, 4), (5, 5), (5, 6), (6, 3), (6, 5), (6, 6)}
}

# Define possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Function to check if a position is within bounds and not a wall
def is_valid(pos, walls):
    return pos not in walls

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes']
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_valid(new_box, state['walls']) and new_box not in boxes:
            new_boxes = boxes - {new_player} | {new_box}
            return {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'walls': state['walls']}
    elif is_valid(new_player, state['walls']):
        return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}
    return None

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state:
                state_id = (new_state['player'], frozenset(new_state['boxes']))
                if state_id not in visited:
                    visited.add(state_id)
                    queue.append((new_state, path + move))
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)