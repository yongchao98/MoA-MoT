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

def is_goal_state(stacks):
    goal = [['A', 'C'], ['B', 'D', 'E', 'G'], ['F'], []]
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i in range(len(stacks)):
        if stacks[i]:  # if stack is not empty
            for j in range(len(stacks)):
                if i != j:  # can't move to same stack
                    moves.append((i+1, j+1))
    return moves

def solve_blocksworld():
    # Initial state
    initial_stacks = [['G'], ['C', 'B'], ['E'], ['A', 'F', 'D']]
    visited_states = set()
    queue = [(initial_stacks, [])]
    
    while queue:
        current_stacks, moves = queue.pop(0)
        
        # Convert current state to tuple for hashing
        state_tuple = tuple(tuple(stack) for stack in current_stacks)
        if state_tuple in visited_states:
            continue
        visited_states.add(state_tuple)
        
        if is_goal_state(current_stacks):
            return moves
        
        possible_moves = get_possible_moves(current_stacks)
        for from_stack, to_stack in possible_moves:
            if is_valid_move(current_stacks, from_stack, to_stack):
                new_stacks = deepcopy(current_stacks)
                block = make_move(new_stacks, from_stack, to_stack)
                new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                queue.append((new_stacks, new_moves))
    
    return None

# Solve and print the solution
solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))