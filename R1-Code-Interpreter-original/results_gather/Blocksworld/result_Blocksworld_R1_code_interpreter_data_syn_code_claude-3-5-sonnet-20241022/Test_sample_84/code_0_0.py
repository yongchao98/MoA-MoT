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

def is_goal_state(stacks, goal):
    return stacks == goal

def find_solution():
    # Initial state
    initial = [['E', 'C'], ['D'], ['A', 'B']]
    goal = [['B'], ['A', 'C', 'E'], ['D']]
    
    queue = deque([(initial, [])])
    visited = {state_to_string(initial)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return path
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    if is_valid_move(current_state, from_stack, to_stack):
                        new_state = make_move(current_state, from_stack, to_stack)
                        state_str = state_to_string(new_state)
                        
                        if state_str not in visited:
                            visited.add(state_str)
                            block = current_state[from_stack-1][-1]
                            move = f"Move {block} from {from_stack} to {to_stack}"
                            queue.append((new_state, path + [move]))

solution = find_solution()
print("<<<" + "\n".join(solution) + ">>>")