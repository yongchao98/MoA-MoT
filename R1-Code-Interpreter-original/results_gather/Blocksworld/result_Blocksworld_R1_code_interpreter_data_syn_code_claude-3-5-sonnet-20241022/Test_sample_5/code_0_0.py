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
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return stacks

def get_state_string(stacks):
    return str(stacks)

def is_goal_state(current, goal):
    return current == goal

def find_solution():
    # Initial state
    initial_stacks = [['D', 'C', 'F', 'B'], ['E'], [], ['G', 'A']]
    goal_stacks = [['B', 'F', 'G'], ['C', 'E'], ['A', 'D'], []]
    
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
            
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack:
                    if is_valid_move(copy.deepcopy(current_stacks), from_stack, to_stack):
                        new_stacks = copy.deepcopy(current_stacks)
                        block = new_stacks[from_stack-1][-1]
                        make_move(new_stacks, from_stack, to_stack)
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))
    
    return None

solution = find_solution()
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")