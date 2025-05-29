import json

def get_next_move(current_pos, target_pos, occupied):
    # Convert positions to coordinates
    cx, cy = map(int, current_pos.replace('C', '').split(','))
    tx, ty = map(int, target_pos.replace('C', '').split(','))
    
    # Try moving in the direction of target
    possible_moves = []
    if cx < tx:
        possible_moves.append(f"C{cx+1},{cy}")
    elif cx > tx:
        possible_moves.append(f"C{cx-1},{cy}")
    if cy < ty:
        possible_moves.append(f"C{cx},{cy+1}")
    elif cy > ty:
        possible_moves.append(f"C{cx},{cy-1}")
    
    # Filter valid moves
    for move in possible_moves:
        if move not in occupied and 1 <= int(move.split(',')[0][1]) <= 5 and 1 <= int(move.split(',')[1]) <= 5:
            return move
    return None

def solve_puzzle():
    initial_state = {
        "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
        "box4": "C4,5", "box5": "C4,4"
    }
    goal_state = {
        "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
        "box4": "C3,1", "box5": "C1,3"
    }
    
    current_state = initial_state.copy()
    solution = [current_state.copy()]
    
    # Move sequence to avoid conflicts
    move_sequence = ["box5", "box2", "box4", "box3", "box1"]
    
    for box in move_sequence:
        while current_state[box] != goal_state[box]:
            # Get occupied positions excluding current box
            occupied = {pos for b, pos in current_state.items() if b != box}
            
            # Get next move
            next_pos = get_next_move(current_state[box], goal_state[box], occupied)
            
            if next_pos:
                current_state[box] = next_pos
                solution.append(current_state.copy())
            else:
                # If stuck, try alternative path
                current_pos = current_state[box]
                cx, cy = map(int, current_pos.replace('C', '').split(','))
                
                # Try moving to any available adjacent cell
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    new_x, new_y = cx + dx, cy + dy
                    if 1 <= new_x <= 5 and 1 <= new_y <= 5:
                        next_pos = f"C{new_x},{new_y}"
                        if next_pos not in occupied:
                            current_state[box] = next_pos
                            solution.append(current_state.copy())
                            break
    
    return solution

# Generate and print solution
solution = solve_puzzle()
print(json.dumps(solution))