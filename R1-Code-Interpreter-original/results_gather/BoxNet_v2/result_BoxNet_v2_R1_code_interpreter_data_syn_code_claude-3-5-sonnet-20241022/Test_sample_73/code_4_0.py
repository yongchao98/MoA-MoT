import json
from collections import deque

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

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    # Sort boxes by distance to goal
    boxes = sorted(state.keys(), 
                  key=lambda b: manhattan_distance(state[b], goal_state[b]),
                  reverse=True)
    
    # Try to move boxes that are not in their goal positions
    for box in boxes:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        # Sort adjacent positions by distance to goal
        next_positions = sorted(adjacency[current_pos],
                              key=lambda pos: manhattan_distance(pos, goal_state[box]))
        
        for next_pos in next_positions:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def solve():
    visited = set()
    queue = [(initial_state, [initial_state])]
    
    while queue:
        current_state, path = queue.pop(0)
        state_str = str(current_state)
        
        if current_state == goal_state:
            return path
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        # Get valid moves and sort them by total distance to goals
        next_states = get_valid_moves(current_state)
        next_states.sort(key=lambda s: sum(manhattan_distance(s[b], goal_state[b]) for b in s))
        
        # Add new states to queue
        for next_state in next_states:
            if str(next_state) not in visited:
                new_path = path + [next_state]
                queue.append((next_state, new_path))
        
        # Keep queue size manageable
        if len(queue) > 1000:
            queue = queue[:500]
        
        # Sort queue by total distance to goals
        queue.sort(key=lambda x: sum(manhattan_distance(x[0][b], goal_state[b]) for b in x[0]))
    
    return None

# Find and print solution
solution = solve()
if solution:
    # Only keep states that represent actual changes
    filtered_solution = [solution[0]]
    for i in range(1, len(solution)):
        if solution[i] != solution[i-1]:
            filtered_solution.append(solution[i])
    print(json.dumps(filtered_solution))
else:
    print("No solution found")