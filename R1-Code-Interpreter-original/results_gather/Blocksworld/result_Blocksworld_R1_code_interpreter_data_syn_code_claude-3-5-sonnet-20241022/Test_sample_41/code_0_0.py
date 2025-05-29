from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty and destination exists
    return (len(stacks[from_stack-1]) > 0 and 
            to_stack <= len(stacks))

def get_top_block(stacks, stack_num):
    if len(stacks[stack_num-1]) > 0:
        return stacks[stack_num-1][-1]
    return None

def make_move(stacks, from_stack, to_stack):
    if is_valid_move(stacks, from_stack, to_stack):
        block = stacks[from_stack-1].pop()
        stacks[to_stack-1].append(block)
        return True
    return False

def state_to_string(stacks):
    return ';'.join([','.join(stack) for stack in stacks])

def is_goal_state(current, goal):
    return state_to_string(current) == state_to_string(goal)

def find_solution():
    # Initial state
    initial_state = [
        ['I'],
        ['G', 'C', 'H', 'A', 'E'],
        ['F', 'B'],
        ['D', 'K', 'J']
    ]
    
    # Goal state
    goal_state = [
        ['A', 'D', 'F', 'I', 'J'],
        ['G', 'H'],
        ['B', 'C', 'E', 'K'],
        []
    ]
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        current_state_str = state_to_string(current_state)
        
        if current_state_str in visited:
            continue
            
        visited.add(current_state_str)
        
        if is_goal_state(current_state, goal_state):
            return moves
            
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack:
                    if is_valid_move(current_state, from_stack, to_stack):
                        new_state = copy.deepcopy(current_state)
                        block = get_top_block(current_state, from_stack)
                        if make_move(new_state, from_stack, to_stack):
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_state, new_moves))

solution = find_solution()
print('\n'.join(solution))