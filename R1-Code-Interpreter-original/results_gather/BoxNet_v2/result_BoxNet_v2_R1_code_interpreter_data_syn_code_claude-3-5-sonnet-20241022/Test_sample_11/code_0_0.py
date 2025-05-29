import json
from heapq import heappush, heappop
from copy import deepcopy

# Problem setup
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
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_heuristic(state):
    h = 0
    for box in state:
        if state[box] != goal_state[box]:
            h += manhattan_distance(state[box], goal_state[box])
    return h

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_valid_moves(state):
    moves = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            # Check if next_pos is occupied by another box
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def solve():
    start_state = initial_state.copy()
    visited = set()
    queue = [(get_heuristic(start_state), 0, [start_state], start_state)]
    
    while queue:
        _, cost, path, current_state = heappop(queue)
        
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
                heappush(queue, (new_cost + get_heuristic(next_state), new_cost, new_path, next_state))
    
    return None

# Solve and print the solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")