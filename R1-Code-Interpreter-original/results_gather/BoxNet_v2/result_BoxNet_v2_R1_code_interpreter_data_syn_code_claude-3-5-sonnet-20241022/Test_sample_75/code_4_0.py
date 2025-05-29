from queue import PriorityQueue
import json
import copy

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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def find_solution():
    # Helper function to check if move is valid
    def is_valid_move(state, box, new_pos):
        occupied = set(state.values())
        occupied.remove(state[box])
        return new_pos not in occupied

    # Helper function to get direct path
    def get_path(start_pos, end_pos, blocked_positions):
        if start_pos == end_pos:
            return []
        
        visited = {start_pos}
        queue = [(start_pos, [])]
        
        while queue:
            pos, path = queue.pop(0)
            for next_pos in adjacency[pos]:
                if next_pos not in visited and next_pos not in blocked_positions:
                    if next_pos == end_pos:
                        return path + [next_pos]
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    # Initialize solution path
    solution = [initial_state]
    current_state = copy.deepcopy(initial_state)
    
    # Move boxes one by one to their goals
    boxes_order = ["box1", "box2", "box3", "box4"]
    
    for box in boxes_order:
        if current_state[box] == goal_state[box]:
            continue
            
        # Get other boxes' positions
        blocked = set(current_state[b] for b in current_state if b != box)
        
        # Find path to goal
        path = get_path(current_state[box], goal_state[box], blocked)
        if path:
            # Move box along path
            for new_pos in path:
                new_state = copy.deepcopy(current_state)
                new_state[box] = new_pos
                current_state = new_state
                solution.append(current_state)
    
    return solution

# Generate and print solution
solution = find_solution()
print(json.dumps(solution))