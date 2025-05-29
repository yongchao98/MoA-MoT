import json
from heapq import heappush, heappop

initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

class PQItem:
    def __init__(self, priority, step, state, path):
        self.priority = priority
        self.step = step
        self.state = state
        self.path = path
    
    def __lt__(self, other):
        return (self.priority, self.step) < (other.priority, other.step)

def get_adjacent(pos):
    r, c = int(pos[1]), int(pos[3])
    adj = []
    if r > 1: adj.append(f"C{r-1},{c}")
    if r < 4: adj.append(f"C{r+1},{c}")
    if c > 1: adj.append(f"C{r},{c-1}")
    if c < 8: adj.append(f"C{r},{c+1}")
    return adj

def manhattan(pos1, pos2):
    r1, c1 = int(pos1[1]), int(pos1[3])
    r2, c2 = int(pos2[1]), int(pos2[3])
    return abs(r1 - r2) + abs(c1 - c2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_next_states(state):
    states = []
    # Sort boxes by distance to goal to prioritize moving boxes that are far from goals
    boxes = sorted(state.keys(), 
                  key=lambda b: manhattan(state[b], goal_state[b]) if state[b] != goal_state[b] else -1,
                  reverse=True)
    
    for box in boxes:
        if state[box] == goal_state[box]:
            continue
        
        pos = state[box]
        for new_pos in get_adjacent(pos):
            if new_pos not in state.values():
                next_state = state.copy()
                next_state[box] = new_pos
                states.append(next_state)
    return states

def solve():
    visited = set()
    queue = []
    
    # Initial state
    h_score = sum(manhattan(initial_state[b], goal_state[b]) for b in initial_state)
    heappush(queue, PQItem(h_score, 0, initial_state, [initial_state]))
    visited.add(state_to_string(initial_state))
    
    while queue:
        current = heappop(queue)
        
        if current.state == goal_state:
            return current.path
            
        if current.step > 20:  # Limit path length
            continue
            
        for next_state in get_next_states(current.state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h_score = sum(manhattan(next_state[b], goal_state[b]) for b in next_state)
                new_path = current.path + [next_state]
                heappush(queue, PQItem(current.step + 1 + h_score, 
                                     current.step + 1,
                                     next_state,
                                     new_path))
    return None

# Find solution
solution = solve()
if solution:
    # Format solution
    result = []
    for state in solution:
        result.append(dict(sorted(state.items())))
    print(json.dumps(result))
else:
    print("No solution found")