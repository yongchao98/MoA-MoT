from collections import deque
import copy

def get_top_blocks(state):
    top_blocks = []
    for i, stack in enumerate(state):
        if stack:
            top_blocks.append((stack[-1], i))
    return top_blocks

def is_valid_move(from_stack, to_stack, state):
    # Check if from_stack has blocks
    if not state[from_stack]:
        return False
    return True

def make_move(state, from_stack, to_stack):
    if not is_valid_move(from_stack, to_stack, state):
        return None
    new_state = copy.deepcopy(state)
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def state_to_string(state):
    return str(state)

def solve_blocksworld():
    # Initial state
    initial_state = [['F'], ['B', 'D'], ['C', 'E', 'A'], []]
    # Goal state
    goal_state = [['D', 'E'], ['B', 'C'], ['A', 'F'], []]
    
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if current_state == goal_state:
            return moves
            
        # Get all possible moves
        for block, from_stack in get_top_blocks(current_state):
            for to_stack in range(4):
                if from_stack != to_stack:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        new_moves = moves + [f"Move {block} from {from_stack+1} to {to_stack+1}"]
                        queue.append((new_state, new_moves))
    
    return None

# Solve and print the solution
solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")