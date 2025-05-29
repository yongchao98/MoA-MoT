import heapq
import json
from copy import deepcopy

# Initial and goal states
initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Simplified adjacency
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
    def __init__(self, priority, moves, state, path):
        self.priority = priority
        self.moves = moves
        self.state = state
        self.path = path

    def __lt__(self, other):
        return (self.priority, self.moves) < (other.priority, other.moves)

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) 
              for box in state if state[box] != goal_state[box])

def is_valid_move(current_pos, new_pos, current_state):
    if new_pos not in adjacency[current_pos]:
        return False
    return new_pos not in current_state.values()

def solve():
    visited = set()
    start_h = get_heuristic(initial_state)
    pq = [PrioritizedState(start_h, 0, initial_state, [initial_state])]
    heapq.heapify(pq)
    
    while pq:
        current = heapq.heappop(pq)
        current_state = current.state
        moves = current.moves
        path = current.path
        
        if current_state == goal_state:
            return path
            
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
        visited.add(state_str)
        
        # Try moving each box that's not at its goal
        for box in current_state:
            if current_state[box] == goal_state[box]:
                continue
                
            current_pos = current_state[box]
            for new_pos in adjacency[current_pos]:
                if is_valid_move(current_pos, new_pos, current_state):
                    new_state = deepcopy(current_state)
                    new_state[box] = new_pos
                    
                    if state_to_string(new_state) not in visited:
                        h = get_heuristic(new_state)
                        new_item = PrioritizedState(
                            h + moves + 1,
                            moves + 1,
                            new_state,
                            path + [new_state]
                        )
                        heapq.heappush(pq, new_item)
    
    return None

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")