import json
from collections import deque
import copy

# Problem setup
initial_state = {
    "box1": "C1,3",
    "box2": "C2,4",
    "box3": "C1,6",
    "box4": "C2,2",
    "box5": "C1,4"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,2",
    "box3": "C1,5",
    "box4": "C2,2",
    "box5": "C2,6"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C1,6"]
}

def state_to_string(state):
    return str(sorted(state.items()))

def is_valid_move(current_state, box, new_pos):
    # Check if the new position is occupied
    for other_box, pos in current_state.items():
        if other_box != box and pos == new_pos:
            return False
    return True

def get_manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_total_distance(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:  # Only count distance if not at goal
            total += get_manhattan_distance(pos, goal_state[box])
    return total

def find_solution():
    start_state = initial_state
    queue = deque([(start_state, [start_state])])
    visited = {state_to_string(start_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we reached the goal
        if current_state == goal_state:
            return path
            
        # Try moving each box
        for box in current_state:
            # Skip if box is already at goal
            if current_state[box] == goal_state[box]:
                continue
                
            current_pos = current_state[box]
            # Try each adjacent position
            for new_pos in adjacency[current_pos]:
                if is_valid_move(current_state, box, new_pos):
                    new_state = copy.deepcopy(current_state)
                    new_state[box] = new_pos
                    
                    state_str = state_to_string(new_state)
                    if state_str not in visited:
                        visited.add(state_str)
                        queue.append((new_state, path + [new_state]))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")