import json
from heapq import heappush, heappop
from copy import deepcopy

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

class StateNode:
    def __init__(self, state, path, g_score):
        self.state = state
        self.path = path
        self.g_score = g_score
        self.h_score = sum(manhattan_distance(state[box], goal_state[box]) 
                          for box in state if state[box] != goal_state[box])
        self.f_score = self.g_score + self.h_score
    
    def __lt__(self, other):
        return (self.f_score, self.g_score) < (other.f_score, other.g_score)

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

def get_next_states(state):
    states = []
    occupied = set(state.values())
    
    # Try to move boxes that are not in their goal positions first
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        for next_pos in adjacency[state[box]]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                states.append(new_state)
    
    return states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def find_solution():
    visited = set()
    queue = [StateNode(initial_state, [initial_state], 0)]
    
    while queue:
        current = heappop(queue)
        
        if current.state == goal_state:
            return current.path
            
        state_str = state_to_string(current.state)
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_next_states(current.state):
            if state_to_string(next_state) not in visited:
                new_path = current.path + [next_state]
                heappush(queue, StateNode(next_state, new_path, current.g_score + 1))
        
        # Break if path is too long
        if len(visited) > 1000:
            break
    
    return None

# Find solution and format output
solution = find_solution()
if solution:
    # Verify solution
    if solution[-1] == goal_state:
        # Format with <<<>>> markers
        result = "<<<" + json.dumps(solution) + ">>>"
        print(result)
    else:
        print("No complete solution found")
else:
    print("No solution found")