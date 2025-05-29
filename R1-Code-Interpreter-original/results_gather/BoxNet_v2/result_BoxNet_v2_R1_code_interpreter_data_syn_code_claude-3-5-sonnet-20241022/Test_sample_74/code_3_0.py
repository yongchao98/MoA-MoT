import heapq
import json

INITIAL_STATE = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}
GOAL_STATE = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

ADJACENCY = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C3,6"]
}

def manhattan_distance(pos1, pos2):
    r1, c1 = map(int, pos1[1:].split(','))
    r2, c2 = map(int, pos2[1:].split(','))
    return abs(r1 - r2) + abs(c1 - c2)

def get_valid_moves(state):
    valid_moves = []
    occupied_cells = set(state.values())
    
    for box, pos in state.items():
        if pos == GOAL_STATE[box]:
            continue
            
        for next_pos in ADJACENCY[pos]:
            if next_pos not in occupied_cells:
                new_state = state.copy()
                new_state[box] = next_pos
                valid_moves.append(new_state)
    
    return valid_moves

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def heuristic(state):
    return sum(manhattan_distance(state[box], GOAL_STATE[box]) for box in state)

def find_solution_greedy():
    current_state = INITIAL_STATE.copy()
    path = [current_state.copy()]
    
    while current_state != GOAL_STATE:
        best_move = None
        best_score = float('inf')
        
        for box in current_state:
            if current_state[box] == GOAL_STATE[box]:
                continue
                
            current_pos = current_state[box]
            for next_pos in ADJACENCY[current_pos]:
                if next_pos not in current_state.values():
                    new_state = current_state.copy()
                    new_state[box] = next_pos
                    score = heuristic(new_state)
                    
                    if score < best_score:
                        best_score = score
                        best_move = (box, next_pos)
        
        if best_move:
            box, next_pos = best_move
            current_state[box] = next_pos
            path.append(current_state.copy())
        else:
            return None
            
    return path

# Find solution and print each state individually
solution = find_solution_greedy()
if solution:
    states_list = []
    for state in solution:
        states_list.append(dict(sorted(state.items())))
    print(json.dumps(states_list))