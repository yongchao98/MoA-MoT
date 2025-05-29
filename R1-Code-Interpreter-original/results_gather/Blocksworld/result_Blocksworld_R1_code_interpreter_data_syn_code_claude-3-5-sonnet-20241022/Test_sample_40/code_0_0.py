from copy import deepcopy
from collections import deque

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def get_state_string(stacks):
    return ';'.join([''.join(stack) for stack in stacks])

def find_solution():
    # Initial state
    initial_stacks = [['G', 'I', 'D'], ['E', 'A', 'K'], ['C', 'B', 'H', 'J', 'F']]
    goal_stacks = [['A', 'D', 'F'], ['B', 'H', 'I', 'J', 'K'], ['C', 'E', 'G']]
    
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = get_state_string(current_stacks)
        
        if current_state == get_state_string(goal_stacks):
            return moves
            
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack and is_valid_move(deepcopy(current_stacks), from_stack, to_stack):
                    new_stacks = deepcopy(current_stacks)
                    block = make_move(new_stacks, from_stack, to_stack)
                    new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                    queue.append((new_stacks, new_moves))

# Find and print solution
solution = find_solution()
print('\n'.join(solution))