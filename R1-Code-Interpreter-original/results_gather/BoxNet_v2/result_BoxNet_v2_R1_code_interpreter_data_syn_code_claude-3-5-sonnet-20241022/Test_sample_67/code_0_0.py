import copy
from collections import deque

def get_path(start, goal, adjacency, occupied_cells):
    if start == goal:
        return []
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_cell in adjacency[current]:
            if next_cell not in visited and next_cell not in occupied_cells:
                if next_cell == goal:
                    return path + [next_cell]
                visited.add(next_cell)
                queue.append((next_cell, path + [next_cell]))
    return None

def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state]
    current_state = copy.deepcopy(initial_state)
    boxes = list(initial_state.keys())
    
    # Keep track of boxes that reached their goals
    completed_boxes = set()
    
    while current_state != goal_state:
        moved = False
        
        # Try to move each box that hasn't reached its goal
        for box in boxes:
            if box in completed_boxes:
                continue
                
            if current_state[box] == goal_state[box]:
                completed_boxes.add(box)
                continue
            
            # Get occupied cells excluding current box
            occupied_cells = {pos for b, pos in current_state.items() if b != box}
            
            # Find path for current box
            path = get_path(current_state[box], goal_state[box], adjacency, occupied_cells)
            
            if path and len(path) > 1:
                # Move box one step along the path
                new_state = copy.deepcopy(current_state)
                new_state[box] = path[1]
                plan.append(new_state)
                current_state = new_state
                moved = True
                break
        
        if not moved:
            return None  # No valid solution found
    
    return plan

# Initial and goal states
initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

# Adjacency dictionary (as provided)
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

# Generate and print the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)