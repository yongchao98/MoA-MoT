import json
from copy import deepcopy
from collections import deque

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

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        # Skip if box is at goal
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def solve():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if current state is goal state
        if current_state == goal_state:
            return path
        
        # Get all possible moves from current state
        for next_state in get_valid_moves(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                
                # Check if this move makes progress towards goal
                progress = False
                for box in next_state:
                    if next_state[box] == goal_state[box] and current_state[box] != goal_state[box]:
                        progress = True
                        break
                
                # Prioritize moves that make progress
                if progress:
                    queue.appendleft((next_state, new_path))
                else:
                    queue.append((next_state, new_path))
                
                # Limit search space
                if len(visited) > 10000:
                    return None
    
    return None

def simplify_solution(solution):
    if not solution:
        return None
    
    # Keep only states where actual changes happen
    simplified = [solution[0]]
    for i in range(1, len(solution)):
        if solution[i] != solution[i-1]:
            simplified.append(solution[i])
    
    return simplified

# Solve and print result
solution = solve()
if solution:
    simplified = simplify_solution(solution)
    print("<<<" + json.dumps(simplified) + ">>>")
else:
    print("No solution found")