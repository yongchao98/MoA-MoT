import heapq
from copy import deepcopy
import json

# Initial and goal states
initial_state = {
    "box1": "C1,2",
    "box2": "C2,4",
    "box3": "C2,3",
    "box4": "C1,4"
}

goal_state = {
    "box1": "C3,1",
    "box2": "C2,3",
    "box3": "C3,2",
    "box4": "C1,5"
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
    "C3,1": ["C3,2", "C2,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5"],
    "C3,6": ["C3,5", "C2,6"]
}

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def heuristic(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:
            total += manhattan_distance(pos, goal_state[box])
    return total

def get_valid_moves(state):
    moves = []
    for box, pos in state.items():
        if pos == goal_state[box]:
            continue
        
        # Get adjacent cells
        for next_pos in adjacency[pos]:
            # Check if the adjacent cell is occupied
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def find_path():
    start_state = initial_state
    visited = set()
    pq = [(heuristic(start_state), 0, [start_state])]
    heapq.heapify(pq)
    
    while pq:
        _, cost, path = heapq.heappop(pq)
        current_state = path[-1]
        
        if all(current_state[box] == goal_state[box] for box in current_state):
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

# Find and print the solution
solution = find_path()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")