from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def state_to_string(stacks):
    return str(stacks)

def is_goal_state(current, goal):
    return current == goal

def find_solution():
    initial = [['A', 'F', 'D', 'E', 'H'], ['B', 'G'], ['C']]
    goal = [['G'], ['A', 'B', 'D', 'E'], ['C', 'F', 'H']]
    
    queue = deque([(initial, [])])
    visited = {state_to_string(initial)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = state_to_string(new_state)
                        if state_str not in visited:
                            visited.add(state_str)
                            block = current_state[from_stack-1][-1]
                            moves_copy = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_state, moves_copy))
    return None

solution = find_solution()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")