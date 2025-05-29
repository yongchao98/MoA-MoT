from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 2),  # Box on goal
    (3, 5),  # Box on goal
    (4, 4),  # Box
    (4, 5),  # Box
    (4, 6)   # Player
)

# Define the goal positions
goals = {(2, 4), (5, 2)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 1), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7),
    (2, 0), (2, 1), (2, 2), (2, 3), (2, 5), (2, 6), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 7),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal(state):
    # Check if all boxes are on goals
    return all(pos in goals for pos in state[:-1])

def is_valid_move(state, move):
    # Check if the move is valid
    player_pos = state[-1]
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    if new_player_pos in walls:
        return False
    
    if new_player_pos in state[:-1]:
        # If the player is moving into a box, check if the box can be moved
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        if new_box_pos in walls or new_box_pos in state[:-1]:
            return False
    
    return True

def apply_move(state, move):
    # Apply the move to the state
    player_pos = state[-1]
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    new_state = list(state)
    new_state[-1] = new_player_pos
    
    if new_player_pos in state[:-1]:
        # Move the box
        box_index = state.index(new_player_pos)
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        new_state[box_index] = new_box_pos
    
    return tuple(new_state)

def solve_sokoban(initial_state):
    # Use BFS to find the solution
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        for move_name, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_name))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)