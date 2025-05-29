import json

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

def is_occupied(pos, state, exclude_box=None):
    return any(box != exclude_box and state[box] == pos for box in state)

def get_valid_moves(state, box):
    current_pos = state[box]
    valid_moves = []
    for next_pos in adjacency[current_pos]:
        if not is_occupied(next_pos, state, box):
            valid_moves.append(next_pos)
    return valid_moves

def solve():
    solution = [dict(initial_state)]
    current_state = dict(initial_state)
    
    # Priority order for moving boxes
    box_order = ["box5", "box2", "box4", "box3", "box1"]
    
    while current_state != goal_state:
        moved = False
        for box in box_order:
            if current_state[box] == goal_state[box]:
                continue
                
            valid_moves = get_valid_moves(current_state, box)
            for move in valid_moves:
                new_state = dict(current_state)
                new_state[box] = move
                
                # Check if this move brings us closer to the goal
                current_x, current_y = map(int, current_state[box].replace('C', '').split(','))
                goal_x, goal_y = map(int, goal_state[box].replace('C', '').split(','))
                move_x, move_y = map(int, move.replace('C', '').split(','))
                
                current_dist = abs(current_x - goal_x) + abs(current_y - goal_y)
                new_dist = abs(move_x - goal_x) + abs(move_y - goal_y)
                
                if new_dist < current_dist or (box == "box5" and move == goal_state[box]):
                    current_state = new_state
                    solution.append(new_state)
                    moved = True
                    break
            
            if moved:
                break
        
        if not moved:
            # If no good moves found, try any valid move
            for box in box_order:
                if current_state[box] != goal_state[box]:
                    valid_moves = get_valid_moves(current_state, box)
                    if valid_moves:
                        new_state = dict(current_state)
                        new_state[box] = valid_moves[0]
                        current_state = new_state
                        solution.append(new_state)
                        break
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))