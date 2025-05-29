import json
from heapq import heappush, heappop
from copy import deepcopy

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
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def is_goal_reached(state):
    return all(state[box] == goal_state[box] for box in state)

def solve():
    open_set = [(get_heuristic(initial_state), 0, state_to_string(initial_state), [initial_state])]
    closed_set = set()
    
    while open_set:
        _, g, state_str, path = heappop(open_set)
        current_state = path[-1]
        
        if is_goal_reached(current_state):
            return path
            
        if state_str in closed_set:
            continue
            
        closed_set.add(state_str)
        
        for next_state in get_valid_moves(current_state):
            next_str = state_to_string(next_state)
            if next_str not in closed_set:
                new_g = g + 1
                new_h = get_heuristic(next_state)
                new_path = path + [next_state]
                heappush(open_set, (new_g + new_h, new_g, next_str, new_path))
    
    return None

solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")