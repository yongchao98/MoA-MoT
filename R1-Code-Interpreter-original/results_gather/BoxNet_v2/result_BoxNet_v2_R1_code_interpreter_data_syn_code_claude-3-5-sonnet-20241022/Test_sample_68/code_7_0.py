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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def get_manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_state_priority(state):
    total_distance = 0
    for box in state:
        if state[box] != goal_state[box]:
            total_distance += get_manhattan_distance(state[box], goal_state[box])
    return total_distance

def get_valid_moves(state):
    moves = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                priority = get_state_priority(new_state)
                moves.append((priority, new_state))
    return sorted(moves)  # Sort by priority

def solve():
    initial = initial_state
    queue = [(get_state_priority(initial), initial, [initial])]
    visited = {json.dumps(initial, sort_keys=True)}
    
    while queue:
        _, current, path = queue.pop(0)
        
        if current == goal_state:
            return path
            
        for priority, next_state in get_valid_moves(current):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((priority, next_state, new_path))
        
        queue.sort()  # Keep queue sorted by priority
    
    return None

# Get solution
solution = solve()

# Format and print solution
if solution:
    # Create list of sorted dictionaries
    formatted_states = []
    for state in solution:
        # Sort each state dictionary by key
        sorted_state = dict(sorted(state.items()))
        formatted_states.append(sorted_state)
    
    # Convert to string with minimal formatting
    result = str(formatted_states).replace("'", '"')
    print(result)
else:
    print("No solution found")