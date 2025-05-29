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

class PriorityEntry:
    def __init__(self, priority, cost, state, path):
        self.priority = priority
        self.cost = cost
        self.state = state
        self.path = path

    def __lt__(self, other):
        return (self.priority, self.cost) < (other.priority, other.cost)

def find_solution():
    start_state = INITIAL_STATE
    visited = set()
    initial_entry = PriorityEntry(heuristic(start_state), 0, start_state, [start_state])
    pq = [initial_entry]
    heapq.heapify(pq)
    
    while pq:
        current = heapq.heappop(pq)
        current_str = state_to_string(current.state)
        
        if current_str in visited:
            continue
        
        visited.add(current_str)
        
        if all(current.state[box] == GOAL_STATE[box] for box in current.state):
            return current.path
            
        for next_state in get_valid_moves(current.state):
            next_str = state_to_string(next_state)
            if next_str not in visited:
                new_path = current.path + [next_state]
                new_cost = current.cost + 1
                new_priority = new_cost + heuristic(next_state)
                heapq.heappush(pq, PriorityEntry(new_priority, new_cost, next_state, new_path))
    
    return None

# Find solution and print in smaller chunks
solution = find_solution()
if solution:
    result = []
    for state in solution:
        result.append(dict(sorted(state.items())))
    print(str(result).replace("'", '"'))