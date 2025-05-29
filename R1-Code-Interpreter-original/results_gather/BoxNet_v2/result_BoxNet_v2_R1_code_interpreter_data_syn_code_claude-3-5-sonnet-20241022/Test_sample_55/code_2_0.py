import json
import heapq
from copy import deepcopy

# Initial and goal states
initial_state = {
    "box1": "C5,5",
    "box2": "C1,3",
    "box3": "C3,3",
    "box4": "C4,5",
    "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5",
    "box2": "C3,4",
    "box3": "C1,5",
    "box4": "C3,1",
    "box5": "C1,3"
}

# Adjacency dictionary
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

class PrioritizedState:
    def __init__(self, priority, steps, state, path):
        self.priority = priority
        self.steps = steps
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.priority, self.steps) < (other.priority, other.steps)

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_valid_moves(state):
    moves = []
    occupied_cells = set(state.values())
    
    for box in state:
        # Skip if box is at goal
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells or next_pos == state[box]:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box])
    return total

def solve():
    start_state = initial_state
    visited = set()
    pq = []
    initial_priority = heuristic(start_state)
    heapq.heappush(pq, PrioritizedState(initial_priority, 0, start_state, [start_state]))
    
    while pq:
        current = heapq.heappop(pq)
        current_state = current.state
        state_str = state_to_string(current_state)
        
        if current_state == goal_state:
            return current.path
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_valid_moves(current_state):
            next_str = state_to_string(next_state)
            if next_str not in visited:
                h = heuristic(next_state)
                priority = current.steps + 1 + h
                heapq.heappush(pq, PrioritizedState(
                    priority,
                    current.steps + 1,
                    next_state,
                    current.path + [next_state]
                ))
    
    return None

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")