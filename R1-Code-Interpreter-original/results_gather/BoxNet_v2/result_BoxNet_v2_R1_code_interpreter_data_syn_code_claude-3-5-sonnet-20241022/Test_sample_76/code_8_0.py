import json
from heapq import heappush, heappop
from copy import deepcopy

class StateNode:
    def __init__(self, f_score, g_score, state_str, state, path):
        self.f_score = f_score
        self.g_score = g_score
        self.state_str = state_str
        self.state = state
        self.path = path
    
    def __lt__(self, other):
        return (self.f_score, self.g_score, self.state_str) < (other.f_score, other.g_score, other.state_str)

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

def is_valid_move(current_state, next_state):
    # Check that only one box moved
    differences = 0
    moved_box = None
    for box in current_state:
        if current_state[box] != next_state[box]:
            differences += 1
            moved_box = box
    if differences != 1:
        return False
    
    # Check that the move is to an adjacent cell
    if moved_box:
        return next_state[moved_box] in adjacency[current_state[moved_box]]
    return False

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def find_solution():
    visited = set()
    start_h = heuristic(initial_state)
    start_node = StateNode(start_h, 0, state_to_string(initial_state), initial_state, [initial_state])
    queue = [start_node]
    
    while queue:
        current = heappop(queue)
        
        if current.state == goal_state:
            # Verify the solution path
            path = current.path
            for i in range(len(path) - 1):
                if not is_valid_move(path[i], path[i + 1]):
                    return None
            return path
            
        if current.state_str in visited:
            continue
            
        visited.add(current.state_str)
        
        for next_state in get_valid_moves(current.state):
            next_str = state_to_string(next_state)
            if next_str not in visited:
                g = current.g_score + 1
                h = heuristic(next_state)
                f = g + h
                new_path = current.path + [next_state]
                next_node = StateNode(f, g, next_str, next_state, new_path)
                heappush(queue, next_node)
    
    return None

# Find solution and format output
solution = find_solution()
if solution:
    # Ensure the solution is complete and valid
    if solution[-1] == goal_state:
        result = "<<<" + json.dumps(solution) + ">>>"
        print(result)
    else:
        print("No complete solution found")
else:
    print("No solution found")