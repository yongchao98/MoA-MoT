import heapq
import json
from copy import deepcopy

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"], "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"], "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C4,7": ["C4,6", "C4,8", "C3,7"], "C4,8": ["C4,7", "C3,8"]
}

class State:
    def __init__(self, boxes, g_score, parent=None):
        self.boxes = boxes
        self.g_score = g_score
        self.parent = parent
        self.h_score = self.calculate_heuristic()
        
    def calculate_heuristic(self):
        total = 0
        for box in self.boxes:
            if self.boxes[box] != goal_state[box]:
                row1, col1 = int(self.boxes[box].split(',')[0][1]), int(self.boxes[box].split(',')[1])
                row2, col2 = int(goal_state[box].split(',')[0][1]), int(goal_state[box].split(',')[1])
                total += abs(row1 - row2) + abs(col1 - col2)
        return total
    
    def __lt__(self, other):
        return (self.g_score + self.h_score) < (other.g_score + other.h_score)
    
    def get_state_string(self):
        return json.dumps(self.boxes, sort_keys=True)
    
    def get_path(self):
        path = []
        current = self
        while current:
            path.append(current.boxes)
            current = current.parent
        return list(reversed(path))

def get_next_states(current_state):
    next_states = []
    for box in current_state.boxes:
        if current_state.boxes[box] == goal_state[box]:
            continue
            
        current_pos = current_state.boxes[box]
        occupied = set(current_state.boxes.values())
        occupied.remove(current_pos)
        
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_boxes = current_state.boxes.copy()
                new_boxes[box] = next_pos
                next_states.append(State(new_boxes, current_state.g_score + 1, current_state))
    
    return next_states

def find_solution():
    initial = State(initial_state, 0)
    open_set = [initial]
    closed_set = {initial.get_state_string()}
    
    while open_set:
        current = heapq.heappop(open_set)
        
        if current.boxes == goal_state:
            return current.get_path()
            
        for next_state in get_next_states(current):
            state_str = next_state.get_state_string()
            if state_str not in closed_set:
                closed_set.add(state_str)
                heapq.heappush(open_set, next_state)
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")