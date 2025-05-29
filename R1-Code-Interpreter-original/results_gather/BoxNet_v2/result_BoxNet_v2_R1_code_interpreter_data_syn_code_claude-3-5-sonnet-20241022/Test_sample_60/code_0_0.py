import json
from collections import deque
import heapq

# Grid dimensions and data
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

def get_valid_moves(current_state):
    valid_moves = []
    occupied_cells = set(current_state.values())
    
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = current_state.copy()
                new_state[box] = next_pos
                valid_moves.append(new_state)
    
    return valid_moves

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box])
    return total

def find_solution():
    start_state = initial_state.copy()
    visited = set()
    pq = [(heuristic(start_state), 0, [start_state])]
    heapq.heapify(pq)
    
    while pq:
        _, cost, path = heapq.heappop(pq)
        current_state = path[-1]
        
        if current_state == goal_state:
            return path
            
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
        visited.add(state_str)
        
        for next_state in get_valid_moves(current_state):
            if state_to_string(next_state) not in visited:
                new_path = path + [next_state]
                new_cost = cost + 1
                heapq.heappush(pq, (new_cost + heuristic(next_state), new_cost, new_path))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")