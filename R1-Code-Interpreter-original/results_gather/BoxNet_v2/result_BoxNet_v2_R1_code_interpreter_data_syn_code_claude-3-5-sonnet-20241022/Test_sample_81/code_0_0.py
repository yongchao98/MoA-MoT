import json
from collections import deque
import heapq

# Helper function to calculate Manhattan distance
def manhattan_distance(pos1_str, pos2_str):
    x1, y1 = map(int, pos1_str[1:].split(','))
    x2, y2 = map(int, pos2_str[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

# Initial and goal states
initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def is_valid_move(current_state, box, new_pos):
    # Check if the new position is occupied by another box
    for other_box, pos in current_state.items():
        if other_box != box and pos == new_pos:
            return False
    return True

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        # If box is at goal, don't move it
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if is_valid_move(current_state, box, next_pos):
                new_state = current_state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def heuristic(state):
    total_distance = 0
    for box in state:
        total_distance += manhattan_distance(state[box], goal_state[box])
    return total_distance

def find_solution():
    start_state = initial_state
    visited = set()
    pq = [(heuristic(start_state), 0, [start_state], start_state)]
    heapq.heapify(pq)
    
    while pq:
        _, cost, path, current_state = heapq.heappop(pq)
        current_string = state_to_string(current_state)
        
        if current_state == goal_state:
            return path
            
        if current_string in visited:
            continue
            
        visited.add(current_string)
        
        for next_state in get_next_states(current_state):
            if state_to_string(next_state) not in visited:
                new_path = path + [next_state]
                new_cost = cost + 1
                heapq.heappush(pq, (new_cost + heuristic(next_state), new_cost, new_path, next_state))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")