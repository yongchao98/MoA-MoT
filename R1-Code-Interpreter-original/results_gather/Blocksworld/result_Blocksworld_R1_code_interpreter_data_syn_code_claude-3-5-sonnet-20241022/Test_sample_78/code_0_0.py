from collections import deque
import copy

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
    initial_stacks = [['F', 'I', 'H'], ['G', 'A', 'E', 'C', 'B'], ['D']]
    goal_stacks = [['A', 'C'], ['E', 'F', 'G'], ['B', 'D', 'H', 'I']]
    
    queue = deque([(initial_stacks, [])])
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        # Check if we reached the goal
        if current_stacks == goal_stacks:
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    if is_valid_move(copy.deepcopy(current_stacks), from_stack, to_stack):
                        new_stacks = copy.deepcopy(current_stacks)
                        block = make_move(new_stacks, from_stack, to_stack)
                        
                        state_string = get_state_string(new_stacks)
                        if state_string not in visited:
                            visited.add(state_string)
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_stacks, new_moves))
    
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")