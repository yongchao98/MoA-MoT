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

def get_state_string(stacks):
    return ';'.join([''.join(stack) for stack in stacks])

def is_goal_state(stacks, goal_stacks):
    return all(stack == goal_stack for stack, goal_stack in zip(stacks, goal_stacks))

def find_solution():
    # Initial state
    initial_stacks = [['G', 'E', 'B'], ['A'], ['C', 'D', 'H', 'F']]
    goal_stacks = [['D', 'F', 'G'], ['A', 'C', 'E', 'H'], ['B']]
    
    queue = deque([(initial_stacks, [])])
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        if is_goal_state(current_stacks, goal_stacks):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        state_string = get_state_string(new_stacks)
                        if state_string not in visited:
                            visited.add(state_string)
                            block = current_stacks[from_stack-1][-1]
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_stacks, new_moves))
    
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")