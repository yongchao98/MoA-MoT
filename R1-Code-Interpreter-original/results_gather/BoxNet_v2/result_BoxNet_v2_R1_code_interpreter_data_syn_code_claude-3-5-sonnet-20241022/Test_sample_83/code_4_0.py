import json
from copy import deepcopy
from collections import deque

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

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    # Try to move boxes that are not in their goal positions
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
    queue = deque([(initial_state, [initial_state])])
    visited = {str(sorted(initial_state.items()))}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we reached the goal state
        if current_state == goal_state:
            return path
            
        # Get all possible moves from current state
        next_states = get_valid_moves(current_state)
        
        # Sort moves by Manhattan distance to goals
        next_states.sort(key=lambda state: sum(
            manhattan_distance(state[box], goal_state[box])
            for box in state
        ))
        
        # Add valid moves to queue
        for next_state in next_states:
            state_str = str(sorted(next_state.items()))
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
                
        # Prevent memory overflow by limiting queue size
        if len(queue) > 5000:
            # Keep states with better heuristic values
            queue = deque(sorted(queue, 
                key=lambda x: sum(manhattan_distance(x[0][box], goal_state[box]) for box in x[0]))[:2500])

# Run the solver and print the solution
solution = solve()
if solution:
    # Format solution for output
    formatted_solution = []
    for state in solution:
        formatted_solution.append({
            "box1": state["box1"],
            "box2": state["box2"],
            "box3": state["box3"],
            "box4": state["box4"],
            "box5": state["box5"]
        })
    print(json.dumps(formatted_solution))
else:
    print("No solution found")