import heapq
import json
from copy import deepcopy

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

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return str(sorted(state.items()))

def get_valid_moves(current_state):
    valid_moves = []
    occupied_cells = set(current_state.values())
    
    for box, pos in current_state.items():
        # If box is at goal, skip it
        if pos == goal_state[box]:
            continue
            
        # Try moving to adjacent cells
        for next_pos in adjacency[pos]:
            if next_pos not in occupied_cells:
                new_state = current_state.copy()
                new_state[box] = next_pos
                valid_moves.append(new_state)
    
    return valid_moves

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) 
              for box in state.keys())

def search_solution():
    start_h = heuristic(initial_state)
    pq = [(start_h, 0, initial_state, [initial_state])]
    visited = {state_to_string(initial_state)}
    
    while pq:
        _, cost, current_state, path = heapq.heappop(pq)
        
        if current_state == goal_state:
            return path
            
        for next_state in get_valid_moves(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h = heuristic(next_state)
                new_path = path + [next_state]
                heapq.heappush(pq, (cost + 1 + h, cost + 1, next_state, new_path))
    
    return None

# Find solution and print
solution = search_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")