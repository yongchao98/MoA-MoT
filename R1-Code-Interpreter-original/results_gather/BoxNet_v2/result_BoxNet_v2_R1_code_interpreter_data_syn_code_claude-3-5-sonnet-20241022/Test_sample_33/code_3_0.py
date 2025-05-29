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

class StateNode:
    def __init__(self, priority, state_id, state, path):
        self.priority = priority
        self.state_id = state_id
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.priority, self.state_id) < (other.priority, other.state_id)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box, pos in state.items():
        if pos == goal_state[box]:
            continue
        for adj in adjacency[pos]:
            if adj not in occupied:
                new_state = deepcopy(state)
                new_state[box] = adj
                moves.append(new_state)
    return moves

def is_goal(state):
    return all(state[box] == goal_state[box] for box in state)

def solve():
    state_counter = 0
    queue = [StateNode(0, state_counter, initial_state, [initial_state])]
    visited = {json.dumps(initial_state, sort_keys=True)}
    
    while queue:
        current_node = heapq.heappop(queue)
        current_state = current_node.state
        current_path = current_node.path
        
        if is_goal(current_state):
            return current_path
            
        for next_state in get_valid_moves(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = current_path + [next_state]
                state_counter += 1
                priority = len(new_path)
                heapq.heappush(queue, StateNode(priority, state_counter, next_state, new_path))
    
    return None

solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")