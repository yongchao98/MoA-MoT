import json
from copy import deepcopy
from heapq import heappush, heappop

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

class PrioritizedState:
    def __init__(self, priority, steps, state, path):
        self.priority = priority
        self.steps = steps
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.priority, self.steps) < (other.priority, other.steps)

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

def get_state_hash(state):
    return tuple(sorted(state.items()))

def heuristic(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:
            total += manhattan_distance(pos, goal_state[box])
    return total

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
                moves.append(new_state)
    
    return moves

def find_solution():
    start_hash = get_state_hash(initial_state)
    visited = {start_hash}
    h_score = heuristic(initial_state)
    queue = [PrioritizedState(h_score, 0, initial_state, [initial_state])]
    
    while queue:
        current = heappop(queue)
        current_state = current.state
        
        if all(current_state[box] == goal_state[box] for box in current_state):
            return current.path
            
        for next_state in get_valid_moves(current_state):
            next_hash = get_state_hash(next_state)
            if next_hash not in visited:
                visited.add(next_hash)
                new_path = current.path + [next_state]
                h_score = heuristic(next_state)
                priority = len(new_path) + h_score
                heappush(queue, PrioritizedState(priority, len(new_path), next_state, new_path))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")