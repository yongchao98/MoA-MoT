import json
from collections import deque

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"], "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"], "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C4,7": ["C4,6", "C4,8", "C3,7"], "C4,8": ["C4,7", "C3,8"]
}

def get_valid_moves(state, target_box):
    moves = []
    occupied = set(state.values())
    current_pos = state[target_box]
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            new_state = state.copy()
            new_state[target_box] = next_pos
            moves.append(new_state)
    
    return moves

def move_box_to_target(current_state, box, target_pos, visited_states):
    queue = deque([(current_state, [current_state])])
    visited = {json.dumps(current_state, sort_keys=True)}
    
    while queue:
        state, path = queue.popleft()
        
        if state[box] == target_pos:
            return path
            
        for next_state in get_valid_moves(state, box):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited and state_str not in visited_states:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
                
        if len(visited) > 1000:  # Limit search space
            return None
            
    return None

def solve():
    current_state = initial_state.copy()
    full_path = [current_state]
    visited_states = {json.dumps(current_state, sort_keys=True)}
    
    # Define box movement order based on dependencies
    box_order = ["box1", "box6", "box3", "box4", "box5", "box2"]
    
    for box in box_order:
        target_pos = goal_state[box]
        if current_state[box] != target_pos:
            path = move_box_to_target(current_state, box, target_pos, visited_states)
            
            if path is None:
                return None
                
            # Update current state and add new states to path
            current_state = path[-1]
            full_path.extend(path[1:])
            
            # Add new states to visited
            for state in path:
                visited_states.add(json.dumps(state, sort_keys=True))
    
    return full_path

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")