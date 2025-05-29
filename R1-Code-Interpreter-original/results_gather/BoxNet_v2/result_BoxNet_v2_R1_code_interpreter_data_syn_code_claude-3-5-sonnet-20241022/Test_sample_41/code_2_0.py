import json
import heapq
from collections import defaultdict

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
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

class PriorityState:
    def __init__(self, f_score, depth, state, path):
        self.f_score = f_score
        self.depth = depth
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.f_score, self.depth) < (other.f_score, other.depth)

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) 
              for box in state if state[box] != goal_state[box])

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        for next_pos in adjacency[current_state[box]]:
            if next_pos in current_state.values():
                continue
            new_state = current_state.copy()
            new_state[box] = next_pos
            next_states.append(new_state)
    return next_states

def solve():
    start = initial_state
    visited = {state_to_string(start)}
    pq = [PriorityState(heuristic(start), 0, start, [start])]
    max_depth = 30
    
    while pq:
        current = heapq.heappop(pq)
        
        if current.state == goal_state:
            return current.path
            
        if current.depth >= max_depth:
            continue
            
        for next_state in get_next_states(current.state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h_score = heuristic(next_state)
                f_score = h_score + current.depth + 1
                next_path = current.path + [next_state]
                heapq.heappush(pq, PriorityState(f_score, current.depth + 1, next_state, next_path))
    
    return None

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")