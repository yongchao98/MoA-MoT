import json
from collections import deque

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

def is_valid_move(state, box, new_pos):
    # Check if new position is occupied by another box
    return new_pos not in state.values()

def get_next_states(state):
    next_states = []
    for box in state:
        # Skip if box is already at goal
        if state[box] == goal_state[box]:
            continue
        
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if is_valid_move(state, box, next_pos):
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def solve_sequential():
    # Solve one box at a time
    current_state = initial_state.copy()
    path = [current_state.copy()]
    
    # Order boxes by priority (you might want to adjust this order)
    box_order = ["box1", "box2", "box3", "box4"]
    
    for box in box_order:
        if current_state[box] == goal_state[box]:
            continue
            
        queue = deque([(current_state.copy(), [])])
        visited = {json.dumps(current_state, sort_keys=True)}
        
        while queue:
            state, moves = queue.popleft()
            
            if state[box] == goal_state[box]:
                # Update current state and path
                for move in moves:
                    path.append(move)
                current_state = state
                break
                
            for next_state in get_next_states(state):
                state_str = json.dumps(next_state, sort_keys=True)
                if state_str not in visited:
                    visited.add(state_str)
                    queue.append((next_state, moves + [next_state]))
    
    return path

# Get solution and print
solution = solve_sequential()
print(json.dumps(solution))