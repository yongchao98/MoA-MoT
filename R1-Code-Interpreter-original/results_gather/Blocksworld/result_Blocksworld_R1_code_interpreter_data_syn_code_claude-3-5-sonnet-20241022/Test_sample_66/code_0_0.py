from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, move_from, move_to):
    if not is_valid_move(stacks, move_from, move_to):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[move_from-1].pop()
    new_stacks[move_to-1].append(block)
    return new_stacks

def state_to_string(stacks):
    return str(stacks)

def is_goal_state(current, goal):
    return current == goal

def find_solution():
    # Initial state
    initial_stacks = [['B', 'A'], ['C', 'F', 'G'], ['D', 'I', 'E', 'J', 'H']]
    goal_stacks = [['A', 'E', 'I'], ['B', 'C', 'J'], ['D', 'F', 'G', 'H']]
    
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = state_to_string(current_stacks)
        
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        if is_goal_state(current_stacks, goal_stacks):
            return moves
            
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        block = current_stacks[from_stack-1][-1]
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))
    
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")