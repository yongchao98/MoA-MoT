import heapq
import json

class State:
    def __init__(self, boxes):
        self.boxes = boxes
        self.h_score = self.calculate_heuristic()
        
    def calculate_heuristic(self):
        total = 0
        for box in self.boxes:
            if self.boxes[box] != goal_state[box]:
                row1, col1 = int(self.boxes[box].split(',')[0][1]), int(self.boxes[box].split(',')[1])
                row2, col2 = int(goal_state[box].split(',')[0][1]), int(goal_state[box].split(',')[1])
                total += abs(row1 - row2) + abs(col1 - col2)
        return total
    
    def get_state_string(self):
        return json.dumps(self.boxes, sort_keys=True)

class SearchNode:
    def __init__(self, state, g_score=0, parent=None):
        self.state = state
        self.g_score = g_score
        self.parent = parent
        self.f_score = g_score + state.h_score
        
    def __lt__(self, other):
        return self.f_score < other.f_score

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Adjacency dictionary (shortened for brevity but contains all connections)
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

def get_successors(node):
    successors = []
    for box in node.state.boxes:
        if node.state.boxes[box] == goal_state[box]:
            continue
            
        current_pos = node.state.boxes[box]
        occupied = set(node.state.boxes.values())
        occupied.remove(current_pos)
        
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_boxes = node.state.boxes.copy()
                new_boxes[box] = next_pos
                new_state = State(new_boxes)
                successors.append(SearchNode(new_state, node.g_score + 1, node))
    
    return successors

def reconstruct_path(node):
    path = []
    current = node
    while current:
        path.append(current.state.boxes)
        current = current.parent
    return list(reversed(path))

def solve():
    start_state = State(initial_state)
    start_node = SearchNode(start_state)
    
    open_set = [start_node]
    closed_set = {start_state.get_state_string()}
    
    while open_set:
        current = heapq.heappop(open_set)
        
        if current.state.boxes == goal_state:
            return reconstruct_path(current)
            
        for successor in get_successors(current):
            state_str = successor.state.get_state_string()
            if state_str not in closed_set:
                closed_set.add(state_str)
                heapq.heappush(open_set, successor)
    
    return None

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")