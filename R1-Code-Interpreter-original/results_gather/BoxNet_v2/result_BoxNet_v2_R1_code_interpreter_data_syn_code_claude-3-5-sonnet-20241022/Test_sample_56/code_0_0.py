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

def is_valid_move(current_state, next_state):
    # Only one box should move at a time
    differences = sum(1 for box in current_state if current_state[box] != next_state[box])
    if differences != 1:
        return False
    
    # The move should be to an adjacent cell
    for box in current_state:
        if current_state[box] != next_state[box]:
            if next_state[box] not in adjacency[current_state[box]]:
                return False
            
    # Check if destination cell is occupied
    destinations = [next_state[box] for box in next_state]
    if len(destinations) != len(set(destinations)):
        return False
    
    return True

def get_possible_moves(current_state):
    moves = []
    for box in current_state:
        # If box is at goal, don't move it
        if current_state[box] == goal_state[box]:
            continue
        
        for adj_cell in adjacency[current_state[box]]:
            # Check if adjacent cell is occupied
            if adj_cell not in current_state.values():
                new_state = copy.deepcopy(current_state)
                new_state[box] = adj_cell
                moves.append(new_state)
    return moves

def manhattan_distance(pos1, pos2):
    row1, col1 = map(int, pos1.split(','))
    row2, col2 = map(int, pos2.split(','))
    return abs(row1 - row2) + abs(col1 - col2)

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def find_solution():
    start_state = initial_state
    visited = {state_to_string(start_state)}
    queue = deque([(start_state, [start_state])])
    
    while queue:
        current_state, path = queue.popleft()
        
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
        
        for next_state in get_possible_moves(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited and is_valid_move(current_state, next_state):
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")