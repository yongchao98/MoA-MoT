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

class State:
    def __init__(self, boxes):
        self.boxes = boxes
    
    def __eq__(self, other):
        return self.boxes == other.boxes
    
    def __hash__(self):
        return hash(tuple(sorted(self.boxes.items())))

class QueueItem:
    def __init__(self, priority, step, path, state):
        self.priority = priority
        self.step = step
        self.path = path
        self.state = state
    
    def __lt__(self, other):
        return (self.priority, self.step) < (other.priority, other.step)

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_next_states(current_boxes):
    next_states = []
    occupied = set(current_boxes.values())
    
    for box in current_boxes:
        if current_boxes[box] == goal_state[box]:
            continue
        
        current_pos = current_boxes[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_boxes = current_boxes.copy()
                new_boxes[box] = next_pos
                next_states.append(new_boxes)
    
    return next_states

def solve():
    initial = State(initial_state)
    goal = State(goal_state)
    visited = {hash(initial)}
    
    # Initial heuristic
    h_score = sum(manhattan_distance(initial_state[box], goal_state[box]) 
                 for box in initial_state)
    
    queue = [QueueItem(h_score, 0, [initial_state], initial_state)]
    heapq.heapify(queue)
    
    while queue:
        current = heapq.heappop(queue)
        current_state = State(current.state)
        
        if current.state == goal_state:
            return current.path
        
        for next_boxes in get_next_states(current.state):
            next_state = State(next_boxes)
            if hash(next_state) not in visited:
                visited.add(hash(next_state))
                h_score = sum(manhattan_distance(next_boxes[box], goal_state[box]) 
                            for box in next_boxes)
                priority = current.step + 1 + h_score
                new_item = QueueItem(
                    priority,
                    current.step + 1,
                    current.path + [next_boxes],
                    next_boxes
                )
                heapq.heappush(queue, new_item)
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")