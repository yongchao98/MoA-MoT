from copy import deepcopy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def get_state_string(stacks):
    return str(stacks)

def is_goal_state(current, goal):
    return current == goal

def find_solution(initial_stacks, goal_stacks):
    from collections import deque
    
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
            
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack and len(current_stacks[from_stack-1]) > 0:
                    new_stacks = deepcopy(current_stacks)
                    block = make_move(new_stacks, from_stack, to_stack)
                    if block:
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))
    
    return None

# Initial state
initial = [['D', 'E'], ['G', 'I', 'F'], ['C', 'B', 'H', 'A']]
# Goal state
goal = [['B', 'C', 'E', 'F', 'I'], ['G', 'H'], ['A', 'D']]

solution = find_solution(initial, goal)
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")