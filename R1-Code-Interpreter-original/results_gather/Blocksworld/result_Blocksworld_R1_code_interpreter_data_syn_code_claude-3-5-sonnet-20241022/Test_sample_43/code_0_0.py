from copy import deepcopy

# Define initial and goal states
initial_state = {
    1: ['D', 'B'],
    2: ['E'],
    3: ['C'],
    4: ['A']
}

goal_state = {
    1: ['D', 'E'],
    2: ['C'],
    3: ['A', 'B']
}

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks
    if not state[from_stack]:
        return False
    return True

def make_move(state, from_stack, to_stack):
    if not is_valid_move(state, from_stack, to_stack):
        return None
    new_state = deepcopy(state)
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def state_to_string(state):
    return str(state)

def is_goal_state(current, goal):
    return all(current[k] == goal[k] for k in goal.keys())

def find_solution(initial, goal):
    visited = set()
    queue = [(initial, [])]
    
    while queue:
        current_state, moves = queue.pop(0)
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal_state(current_state, goal):
            return moves
            
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack and is_valid_move(current_state, from_stack, to_stack):
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        block = get_top_block(current_state[from_stack])
                        move = f"Move {block} from {from_stack} to {to_stack}"
                        queue.append((new_state, moves + [move]))

# Find and print solution
solution = find_solution(initial_state, goal_state)
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")