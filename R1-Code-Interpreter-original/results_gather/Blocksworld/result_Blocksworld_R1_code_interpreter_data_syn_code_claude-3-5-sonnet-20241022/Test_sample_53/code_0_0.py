from collections import deque
import copy

def get_top_blocks(state):
    return [stack[-1] if stack else None for stack in state]

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks and destination exists
    if not state[from_stack-1]:
        return False
    if to_stack > len(state):
        return False
    return True

def make_move(state, from_stack, to_stack):
    if not is_valid_move(state, from_stack, to_stack):
        return None
    new_state = copy.deepcopy(state)
    block = new_state[from_stack-1].pop()
    new_state[to_stack-1].append(block)
    return new_state

def state_to_string(state):
    return ';'.join(','.join(stack) for stack in state)

def find_solution():
    initial_state = [['C'], ['B', 'E', 'A', 'J', 'F'], ['G', 'I'], ['D', 'H']]
    goal_state = [['E'], ['A', 'B', 'C', 'D', 'H'], ['F', 'G', 'I', 'J']]
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = state_to_string(current_state)
        
        if state_str == state_to_string(goal_state):
            return moves
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack:
                    if is_valid_move(current_state, from_stack, to_stack):
                        new_state = make_move(current_state, from_stack, to_stack)
                        if new_state:
                            block = current_state[from_stack-1][-1]
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_state, new_moves))

solution = find_solution()
print('\n'.join(solution))