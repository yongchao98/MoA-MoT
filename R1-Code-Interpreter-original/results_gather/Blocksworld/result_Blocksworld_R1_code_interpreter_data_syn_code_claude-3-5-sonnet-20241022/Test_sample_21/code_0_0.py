from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    if not stacks[from_stack]:
        return False
    return True

def make_move(stacks, move_from, move_to):
    if not is_valid_move(stacks, move_from, move_to):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[move_from].pop()
    new_stacks[move_to].append(block)
    return new_stacks

def state_to_string(stacks):
    return str(stacks)

def is_goal_state(current, goal):
    return current == goal

def find_solution():
    # Initial state
    initial_stacks = [
        [],             # stack1
        ['F'],         # stack2
        ['G', 'C'],    # stack3
        ['A', 'B', 'D', 'E', 'H']  # stack4
    ]
    
    # Goal state
    goal_stacks = [
        ['B', 'D', 'E'],  # stack1
        ['G', 'H'],       # stack2
        ['A', 'C', 'F']   # stack3
    ]
    
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = state_to_string(current_stacks)
        
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        if is_goal_state(current_stacks, goal_stacks):
            return moves
            
        # Try all possible moves
        for from_stack in range(4):
            for to_stack in range(4):
                if from_stack != to_stack and current_stacks[from_stack]:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        block = current_stacks[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_stacks, moves + [move]))

solution = find_solution()
if solution:
    print('\n'.join(solution))