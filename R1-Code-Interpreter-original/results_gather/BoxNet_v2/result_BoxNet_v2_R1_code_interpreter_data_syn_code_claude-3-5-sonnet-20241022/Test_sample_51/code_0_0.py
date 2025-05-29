import json
from queue import PriorityQueue
import copy

# Problem setup
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
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_heuristic(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:
            total += manhattan_distance(pos, goal_state[box])
    return total

def get_valid_moves(state):
    moves = []
    occupied_cells = set(state.values())
    
    for box, current_pos in state.items():
        if current_pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = copy.deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def find_path():
    start_state = initial_state
    visited = set()
    pq = PriorityQueue()
    pq.put((get_heuristic(start_state), [start_state]))
    
    while not pq.empty():
        _, path = pq.get()
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
                pq.put((len(new_path) + get_heuristic(next_state), new_path))
    
    return None

# Find and print solution
solution = find_path()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")