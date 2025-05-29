import json
from heapq import heappush, heappop
from copy import deepcopy

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"], "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"], "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C4,5"]
}

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box, pos in state.items():
        if pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append((new_state, manhattan_distance(next_pos, goal_state[box])))
    
    return moves

def find_solution():
    visited = set()
    start_h = heuristic(initial_state)
    # Priority queue items: (f_score, moves_count, state_str, state, path)
    queue = [(start_h, 0, state_to_string(initial_state), initial_state, [initial_state])]
    
    while queue:
        _, moves, _, current_state, path = heappop(queue)
        
        if current_state == goal_state:
            return path
            
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state, move_cost in get_valid_moves(current_state):
            next_str = state_to_string(next_state)
            if next_str not in visited:
                h = heuristic(next_state)
                g = moves + 1
                f = g + h
                new_path = path + [next_state]
                heappush(queue, (f, g, next_str, next_state, new_path))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")