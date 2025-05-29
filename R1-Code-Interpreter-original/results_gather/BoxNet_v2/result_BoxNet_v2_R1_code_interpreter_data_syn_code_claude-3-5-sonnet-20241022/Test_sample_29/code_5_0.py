import heapq
import json

class PriorityEntry:
    def __init__(self, priority, count, path):
        self.priority = priority
        self.count = count
        self.path = path
    
    def __lt__(self, other):
        return (self.priority, self.count) < (other.priority, other.count)

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_valid_moves(state, adjacency, goal_state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def solve():
    initial_state = {
        "box1": "C1,4", "box2": "C2,8",
        "box3": "C2,2", "box4": "C2,1"
    }
    
    goal_state = {
        "box1": "C2,8", "box2": "C2,2",
        "box3": "C1,2", "box4": "C2,5"
    }
    
    adjacency = {
        "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
        "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
        "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
        "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
        "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
        "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
        "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
        "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
    }
    
    def state_to_string(state):
        return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))
    
    def heuristic(state):
        return sum(manhattan_distance(state[box], goal_state[box]) 
                  for box in state if state[box] != goal_state[box])
    
    visited = set()
    count = 0
    queue = []
    initial_entry = PriorityEntry(heuristic(initial_state), count, [initial_state])
    heapq.heappush(queue, initial_entry)
    
    while queue:
        current_entry = heapq.heappop(queue)
        current_state = current_entry.path[-1]
        
        if current_state == goal_state:
            return current_entry.path
        
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_valid_moves(current_state, adjacency, goal_state):
            if state_to_string(next_state) not in visited:
                count += 1
                new_path = current_entry.path + [next_state]
                priority = len(new_path) + heuristic(next_state)
                heapq.heappush(queue, PriorityEntry(priority, count, new_path))
    
    return None

# Get and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")