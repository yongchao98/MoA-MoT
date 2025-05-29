import heapq
from copy import deepcopy
import json

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

class PQEntry:
    def __init__(self, priority, step_count, path, state):
        self.priority = priority
        self.step_count = step_count
        self.path = path
        self.state = state
        
    def __lt__(self, other):
        return (self.priority, self.step_count) < (other.priority, other.step_count)

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_next_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
        
        curr_pos = state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def solve_boxes():
    visited = set()
    start_h = heuristic(initial_state)
    pq = [PQEntry(start_h, 0, [initial_state], initial_state)]
    heapq.heapify(pq)
    
    while pq:
        current = heapq.heappop(pq)
        current_state = current.state
        current_path = current.path
        steps = current.step_count
        
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if current_state == goal_state:
            return current_path
            
        for next_state in get_next_moves(current_state):
            if state_to_string(next_state) not in visited:
                h = heuristic(next_state)
                new_entry = PQEntry(
                    steps + 1 + h,
                    steps + 1,
                    current_path + [next_state],
                    next_state
                )
                heapq.heappush(pq, new_entry)
    
    return None

solution = solve_boxes()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")