import json
import heapq
from copy import deepcopy

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

class State:
    def __init__(self, boxes_dict):
        self.boxes = boxes_dict

    def __eq__(self, other):
        return self.boxes == other.boxes

    def __hash__(self):
        return hash(tuple(sorted(self.boxes.items())))

    def to_dict(self):
        return self.boxes.copy()

class QueueItem:
    def __init__(self, f_score, g_score, state, path):
        self.f_score = f_score
        self.g_score = g_score
        self.state = state
        self.path = path

    def __lt__(self, other):
        return (self.f_score, self.g_score) < (other.f_score, other.g_score)

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_next_states(current_state):
    next_states = []
    occupied = set(current_state.boxes.values())
    
    for box in current_state.boxes:
        if current_state.boxes[box] == goal_state[box]:
            continue
            
        curr_pos = current_state.boxes[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied or next_pos == current_state.boxes[box]:
                new_boxes = current_state.boxes.copy()
                new_boxes[box] = next_pos
                next_states.append(State(new_boxes))
    
    return next_states

def heuristic(state):
    return sum(manhattan_distance(state.boxes[box], goal_state[box]) 
              for box in state.boxes if state.boxes[box] != goal_state[box])

def solve():
    start_state = State(initial_state)
    goal = State(goal_state)
    visited = set()
    
    # Initialize queue with start state
    start_item = QueueItem(heuristic(start_state), 0, start_state, [start_state.to_dict()])
    queue = [start_item]
    heapq.heapify(queue)
    
    while queue:
        current = heapq.heappop(queue)
        
        if current.state.boxes == goal.boxes:
            return current.path
            
        if current.state in visited:
            continue
            
        visited.add(current.state)
        
        for next_state in get_next_states(current.state):
            if next_state not in visited:
                g_score = current.g_score + 1
                f_score = g_score + heuristic(next_state)
                
                next_path = current.path + [next_state.to_dict()]
                next_item = QueueItem(f_score, g_score, next_state, next_path)
                
                heapq.heappush(queue, next_item)
        
        # Limit search depth
        if current.g_score > 30:
            continue
    
    return None

# Find and print solution
solution = solve()
if solution:
    # Format solution for output
    formatted_solution = []
    for state in solution:
        formatted_solution.append(dict(sorted(state.items())))
    print(json.dumps(formatted_solution))
else:
    print("No solution found")