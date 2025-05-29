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
    initial = [['D', 'C', 'E'], ['A'], [], ['B']]
    goal = [['A', 'C', 'D'], ['B'], ['E'], []]
    
    visited = set()
    queue = deque([(initial, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return moves
        
        state_str = state_to_string(current_state)
        if state_str in visited:
            continue
        
        visited.add(state_str)
        
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack:
                    if is_valid_move(current_state, from_stack, to_stack):
                        block = current_state[from_stack-1][-1]
                        new_state = make_move(current_state, from_stack, to_stack)
                        if new_state:
                            move = f"Move {block} from {from_stack} to {to_stack}"
                            queue.append((new_state, moves + [move]))

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")