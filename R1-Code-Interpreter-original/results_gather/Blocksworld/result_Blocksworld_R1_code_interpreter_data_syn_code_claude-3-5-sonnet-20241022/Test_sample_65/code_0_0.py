from copy import deepcopy
from collections import deque

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = deepcopy(stacks)
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def get_state_string(stacks):
    return ';'.join([','.join(stack) for stack in stacks])

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    # Initial state
    initial_stacks = [
        ['B', 'G', 'I'],
        ['F', 'H', 'D'],
        ['E', 'C'],
        ['A']
    ]
    
    # Goal state
    goal_stacks = [
        ['B', 'D', 'E', 'I'],
        ['F', 'G', 'H'],
        ['A', 'C'],
        []
    ]
    
    visited = set()
    queue = deque([(initial_stacks, [])])
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = get_state_string(current_stacks)
        
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        if is_goal_state(current_stacks, goal_stacks):
            return moves
            
        for from_stack in range(len(current_stacks)):
            if not current_stacks[from_stack]:
                continue
                
            for to_stack in range(len(current_stacks)):
                if from_stack == to_stack:
                    continue
                    
                new_stacks = make_move(current_stacks, from_stack, to_stack)
                if new_stacks:
                    block = current_stacks[from_stack][-1]
                    move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                    queue.append((new_stacks, moves + [move]))
    
    return None

solution = find_solution()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")