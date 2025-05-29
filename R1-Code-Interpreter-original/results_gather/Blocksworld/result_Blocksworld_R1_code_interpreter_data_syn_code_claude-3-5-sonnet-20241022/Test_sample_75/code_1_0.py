from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    if from_stack < 0 or from_stack >= len(state) or not state[from_stack]:
        return False
    return True

def make_move(state, from_stack, to_stack):
    if not is_valid_move(state, from_stack, to_stack):
        return None
    new_state = copy.deepcopy(state)
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def state_to_string(state):
    return ';'.join(','.join(stack) for stack in state)

def find_solution():
    initial_state = [['D', 'C'], ['G'], ['A', 'E'], ['F', 'B']]
    goal_state = [['A', 'E'], ['B', 'D'], ['C', 'F', 'G']]
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        current_str = state_to_string(current_state)
        
        if current_str == state_to_string(goal_state):
            return moves
            
        if current_str in visited:
            continue
            
        visited.add(current_str)
        
        for from_stack in range(len(current_state)):
            if not current_state[from_stack]:
                continue
            for to_stack in range(len(current_state)):
                if from_stack != to_stack:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        block = current_state[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_state, moves + [move]))

solution = find_solution()
if solution:
    print('\n'.join(solution))