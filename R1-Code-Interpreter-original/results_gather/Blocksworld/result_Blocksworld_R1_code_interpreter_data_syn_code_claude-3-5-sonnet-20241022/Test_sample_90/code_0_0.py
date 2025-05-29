from copy import deepcopy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    # Initial state
    initial_state = [['D', 'E', 'G'], ['C', 'B', 'A'], ['F', 'H']]
    goal_state = [['B', 'E', 'G'], ['C', 'D', 'H'], ['A', 'F']]
    
    # BFS search
    queue = [(initial_state, [])]
    visited = set()
    
    while queue:
        current_state, moves = queue.pop(0)
        
        # Convert current state to tuple for hashing
        state_tuple = tuple(tuple(stack) for stack in current_state)
        if state_tuple in visited:
            continue
        visited.add(state_tuple)
        
        if is_goal_state(current_state, goal_state):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    if is_valid_move(deepcopy(current_state), from_stack, to_stack):
                        new_state = deepcopy(current_state)
                        block = make_move(new_state, from_stack, to_stack)
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_state, new_moves))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")