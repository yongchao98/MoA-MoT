import json
from heapq import heappush, heappop
from copy import deepcopy

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

class StateNode:
    def __init__(self, h_score, g_score, state, path):
        self.h_score = h_score
        self.g_score = g_score
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.h_score + self.g_score) < (other.h_score + other.g_score)

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append((box, new_state))
    
    return moves

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def search():
    start_state = initial_state
    visited = set()
    queue = []
    
    initial_h = sum(manhattan_distance(start_state[box], goal_state[box]) for box in start_state)
    heappush(queue, StateNode(initial_h, 0, start_state, [start_state]))
    
    while queue:
        current = heappop(queue)
        current_str = state_to_string(current.state)
        
        if current_str in visited:
            continue
            
        visited.add(current_str)
        
        if all(current.state[box] == goal_state[box] for box in current.state):
            return current.path
            
        for box, next_state in get_valid_moves(current.state):
            next_str = state_to_string(next_state)
            if next_str not in visited:
                h_score = sum(manhattan_distance(next_state[b], goal_state[b]) for b in next_state)
                g_score = current.g_score + 1
                new_path = current.path + [next_state]
                heappush(queue, StateNode(h_score, g_score, next_state, new_path))
    
    return None

# Find and print solution
solution = search()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")