from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (5, 8),
    'boxes': {(3, 3), (4, 5), (4, 6), (4, 7), (5, 7)},
    'goals': {(3, 5), (3, 6), (5, 8)},
    'walls': {(0, i) for i in range(10)} | {(6, i) for i in range(10)} | {(i, 0) for i in range(7)} | {(i, 9) for i in range(7)} | {(i, 8) for i in range(7)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Function to check if a position is free (not a wall or a box)
def is_free(position, boxes, walls):
    return position not in boxes and position not in walls

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes']
    walls = state['walls']
    
    # Calculate new player position
    new_player = (player[0] + moves[move][0], player[1] + moves[move][1])
    
    # Check if the new player position is a box
    if new_player in boxes:
        # Calculate new box position
        new_box = (new_player[0] + moves[move][0], new_player[1] + moves[move][1])
        # Check if the new box position is free
        if is_free(new_box, boxes, walls):
            # Move the box
            boxes = boxes - {new_player} | {new_box}
        else:
            return None
    elif not is_free(new_player, boxes, walls):
        return None
    
    return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        
        if is_solved(state):
            return path
        
        for move in moves:
            new_state = apply_move(state, move)
            if new_state and (new_state['player'], frozenset(new_state['boxes'])) not in visited:
                visited.add((new_state['player'], frozenset(new_state['boxes'])))
                queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)