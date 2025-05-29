import json
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

def get_single_box_moves(state, box):
    moves = []
    current_pos = state[box]
    occupied = set(state.values())
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            new_state = state.copy()
            new_state[box] = next_pos
            moves.append(new_state)
    
    return moves

def solve():
    # Initialize path with starting state
    path = [initial_state]
    current_state = initial_state.copy()
    
    # Order boxes by priority (those that need to move)
    boxes_to_move = []
    for box in current_state:
        if current_state[box] != goal_state[box]:
            boxes_to_move.append(box)
    
    # Move each box to its goal
    for box in boxes_to_move:
        while current_state[box] != goal_state[box]:
            # Get current position and target position
            current_pos = current_state[box]
            target_pos = goal_state[box]
            
            # Find next move towards goal
            best_move = None
            min_distance = float('inf')
            
            for next_state in get_single_box_moves(current_state, box):
                new_pos = next_state[box]
                # Calculate Manhattan distance to goal
                curr_row, curr_col = int(new_pos.split(',')[0][1]), int(new_pos.split(',')[1])
                goal_row, goal_col = int(target_pos.split(',')[0][1]), int(target_pos.split(',')[1])
                distance = abs(curr_row - goal_row) + abs(curr_col - goal_col)
                
                if distance < min_distance:
                    min_distance = distance
                    best_move = next_state
            
            if best_move:
                current_state = best_move
                path.append(current_state)
            else:
                # If no direct move is possible, try moving other boxes out of the way
                for other_box in boxes_to_move:
                    if other_box != box and current_state[other_box] != goal_state[other_box]:
                        for move in get_single_box_moves(current_state, other_box):
                            current_state = move
                            path.append(current_state)
                            break
    
    return path

# Solve and print result
solution = solve()
print(json.dumps(solution))