from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks
    if not state[from_stack-1]:
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
    initial_state = [['A', 'E', 'F'], ['B', 'C'], ['D', 'H', 'G'], ['I']]
    goal_state = [['B', 'C', 'D', 'E'], ['A', 'F', 'G', 'H'], ['I'], []]
    
    queue = deque([(initial_state, [])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack and current_state[from_stack-1]:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = state_to_string(new_state)
                        if state_str not in visited:
                            visited.add(state_str)
                            block = current_state[from_stack-1][-1]
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_state, new_moves))
    return None

solution = find_solution()
if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")