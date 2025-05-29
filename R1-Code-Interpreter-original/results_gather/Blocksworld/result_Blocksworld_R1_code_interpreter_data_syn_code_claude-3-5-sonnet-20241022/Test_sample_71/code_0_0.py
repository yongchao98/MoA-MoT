from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks
    if not state[from_stack]:
        return False
    return True

def get_state_string(state):
    return ','.join([''.join(stack) for stack in state])

def find_solution(initial_state, goal_state):
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = get_state_string(current_state)
        
        if state_str == get_state_string(goal_state):
            return moves
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        # Try all possible moves
        for from_stack in range(4):
            if not current_state[from_stack]:
                continue
            
            block = current_state[from_stack][-1]
            
            for to_stack in range(4):
                if from_stack != to_stack:
                    # Create new state
                    new_state = [list(stack) for stack in current_state]
                    new_state[from_stack].pop()
                    new_state[to_stack].append(block)
                    
                    # Record move
                    move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                    
                    queue.append((new_state, moves + [move]))
    
    return None

# Initial state
initial = [list('CE'), [], list('A'), list('DB')]

# Goal state
goal = [list('DE'), list('BC'), list('A'), []]

# Find solution
solution = find_solution(initial, goal)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")