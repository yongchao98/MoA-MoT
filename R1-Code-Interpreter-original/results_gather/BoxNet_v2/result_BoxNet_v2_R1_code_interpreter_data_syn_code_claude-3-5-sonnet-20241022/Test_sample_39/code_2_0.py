import copy
from collections import deque

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

def get_next_moves(current_pos, occupied_cells):
    return [pos for pos in adjacency[current_pos] if pos not in occupied_cells]

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.split(',')[0][1:]), map(int, pos1.split(',')[1])
    x2, y2 = map(int, pos2.split(',')[0][1:]), map(int, pos2.split(',')[1])
    return abs(x1 - x2) + abs(y1 - y2)

def generate_plan():
    states = [initial_state]
    current_state = copy.deepcopy(initial_state)
    
    while current_state != goal_state:
        moved = False
        # Sort boxes by distance to goal to prioritize moves
        boxes = sorted(current_state.keys(), 
                      key=lambda b: manhattan_distance(current_state[b], goal_state[b]))
        
        for box in boxes:
            if current_state[box] == goal_state[box]:
                continue
                
            occupied_cells = {pos for b, pos in current_state.items() if b != box}
            next_moves = get_next_moves(current_state[box], occupied_cells)
            
            # Find best move towards goal
            best_move = None
            best_distance = float('inf')
            for move in next_moves:
                dist = manhattan_distance(move, goal_state[box])
                if dist < best_distance:
                    best_distance = dist
                    best_move = move
            
            if best_move:
                new_state = copy.deepcopy(current_state)
                new_state[box] = best_move
                current_state = new_state
                states.append(new_state)
                moved = True
                break
        
        if not moved:
            # Try alternative paths if direct path is blocked
            for box in boxes:
                if current_state[box] == goal_state[box]:
                    continue
                    
                occupied_cells = {pos for b, pos in current_state.items() if b != box}
                next_moves = get_next_moves(current_state[box], occupied_cells)
                
                if next_moves:
                    new_state = copy.deepcopy(current_state)
                    new_state[box] = next_moves[0]
                    current_state = new_state
                    states.append(new_state)
                    moved = True
                    break
                    
        if not moved:
            return None
            
        if len(states) > 100:  # Prevent infinite loops
            return None
    
    return states

# Generate and print the plan
plan = generate_plan()
if plan:
    print(plan)
else:
    print("No solution found")