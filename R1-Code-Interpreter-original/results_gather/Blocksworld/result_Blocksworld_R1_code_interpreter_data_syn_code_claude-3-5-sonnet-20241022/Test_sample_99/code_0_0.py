from collections import deque
import copy

def get_top_blocks(state):
    top_blocks = []
    for stack_idx, stack in enumerate(state):
        if stack:  # if stack is not empty
            top_blocks.append((stack[-1], stack_idx))
    return top_blocks

def is_valid_move(block, from_stack, to_stack, state):
    # Check if block is at the top of from_stack
    if not state[from_stack] or state[from_stack][-1] != block:
        return False
    return True

def make_move(state, block, from_stack, to_stack):
    new_state = copy.deepcopy(state)
    new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def state_to_string(state):
    return str(state)

def find_solution():
    initial_state = [['D'], ['E', 'A'], ['C', 'B'], []]
    goal_state = [['A', 'C'], ['B', 'D', 'E'], [], []]
    
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
            
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        # Get all possible moves
        top_blocks = get_top_blocks(current_state)
        for block, from_stack in top_blocks:
            for to_stack in range(4):
                if to_stack != from_stack:
                    if is_valid_move(block, from_stack, to_stack, current_state):
                        new_state = make_move(current_state, block, from_stack, to_stack)
                        new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                        queue.append((new_state, new_moves))
    
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")